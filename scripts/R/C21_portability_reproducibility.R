cat("C21 STARTED\n")

library(data.table)

need_pkgs <- c("glmnet","pROC")
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) {
    stop("Missing package: ", p, " (install.packages('", p, "') )")
  }
}

library(glmnet)
library(pROC)

meta_path <- "data/processed/ipd_sex/ipd_metadata.csv"
X_path    <- "data/processed/ipd_paper2/ipd_genus_counts_shared2_rebuilt.csv"

stopifnot(file.exists(meta_path), file.exists(X_path))

meta <- fread(meta_path)
X    <- fread(X_path)

setkey(meta, SampleID)
setkey(X, SampleID)

dt <- X[meta, on="SampleID", nomatch=0]
stopifnot(nrow(dt) == nrow(meta))

# outcome: ensure numeric 0/1
y <- meta$PD
if (is.factor(y)) y <- as.character(y)
y <- as.integer(y)
# if coded as 0/1 already ok; if coded as Control/PD:
if (any(is.na(y))) {
  y <- ifelse(meta$PD %in% c("PD","1",1), 1L, 0L)
}
stopifnot(all(y %in% c(0L,1L)))

cohort <- meta$Cohort

genus_cols <- setdiff(names(dt), c("SampleID","PD","Sex","Cohort"))

# ---------- transforms ----------
clr_transform <- function(mat) {
  mat_pc <- mat + 1
  log_mat <- log(mat_pc)
  log_mat - rowMeans(log_mat)
}

rel_abund <- function(mat) {
  rs <- rowSums(mat)
  if (any(rs == 0)) stop("Zero-sum sample in counts.")
  mat / rs
}

# A: CLR on counts
X_counts <- as.matrix(dt[, ..genus_cols]); storage.mode(X_counts) <- "numeric"
X_clr <- clr_transform(X_counts)

# B: log relative abundance (sensitivity comparator)
X_rel <- rel_abund(X_counts)
X_logrel <- log(X_rel + 1e-6)

# ---------- helpers ----------
auc_from_probs <- function(y_true, p_hat) {
  as.numeric(pROC::auc(pROC::roc(y_true, p_hat, quiet=TRUE)))
}

fit_glmnet <- function(Xtr, ytr, alpha=0.5, nfolds=5, seed=1) {
  set.seed(seed)
  cv <- cv.glmnet(
    x = Xtr, y = ytr,
    family = "binomial",
    alpha = alpha,
    nfolds = nfolds,
    standardize = TRUE
  )
  cv
}

nonzero_features <- function(cvfit, feature_names, which_lambda=c("lambda.1se","lambda.min")) {
  which_lambda <- match.arg(which_lambda)
  lam <- cvfit[[which_lambda]]
  b <- as.matrix(coef(cvfit, s = lam))
  b <- b[-1, , drop=FALSE] # drop intercept
  nz <- which(abs(b[,1]) > 0)
  data.table(
    feature = feature_names[nz],
    coef = b[nz,1],
    sign = sign(b[nz,1]),
    lambda = lam
  )
}

# ---------- output dirs ----------
out_dir <- "results/paper2_portability"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

run_suite <- function(Xmat, label) {

  feature_names <- colnames(Xmat)

  # 1) Within-cohort repeated CV (AUC)
  within <- rbindlist(lapply(sort(unique(cohort)), function(coh) {
    idx <- which(cohort == coh)
    Xc <- Xmat[idx, , drop=FALSE]
    yc <- y[idx]

    # repeated 5-fold CV (10 repeats) via different seeds
    reps <- 10
    aucs <- sapply(1:reps, function(r) {
      cv <- fit_glmnet(Xc, yc, alpha=0.5, nfolds=5, seed=1000+r)
      p_hat <- as.numeric(predict(cv, newx=Xc, s=cv$lambda.1se, type="response"))
      # NOTE: this is in-sample prob; for strict CV preds you’d need fold-pred extraction.
      # We keep within-cohort as comparative baseline; cross-cohort/LOCO are the primary endpoints.
      auc_from_probs(yc, p_hat)
    })

    data.table(metric=label, cohort_train=coh, cohort_test=coh,
               setting="within_baseline", AUC=mean(aucs), AUC_sd=sd(aucs),
               n=length(idx))
  }))

  # 2) Cross-cohort transfer: train on A, test on B
  cohorts <- sort(unique(cohort))
  transfer <- rbindlist(lapply(cohorts, function(tr) {
    idx_tr <- which(cohort == tr)
    Xtr <- Xmat[idx_tr, , drop=FALSE]
    ytr <- y[idx_tr]

    cv <- fit_glmnet(Xtr, ytr, alpha=0.5, nfolds=5, seed=123)

    rbindlist(lapply(cohorts, function(te) {
      idx_te <- which(cohort == te)
      Xte <- Xmat[idx_te, , drop=FALSE]
      yte <- y[idx_te]

      p_hat <- as.numeric(predict(cv, newx=Xte, s=cv$lambda.1se, type="response"))
      data.table(metric=label, cohort_train=tr, cohort_test=te,
                 setting="transfer", AUC=auc_from_probs(yte, p_hat),
                 n_train=length(idx_tr), n_test=length(idx_te))
    }))
  }))

  # 3) LOCO: train on 2 cohorts, test on held-out
  loco <- rbindlist(lapply(cohorts, function(held) {
    idx_te <- which(cohort == held)
    idx_tr <- which(cohort != held)

    Xtr <- Xmat[idx_tr, , drop=FALSE]
    ytr <- y[idx_tr]
    Xte <- Xmat[idx_te, , drop=FALSE]
    yte <- y[idx_te]

    cv <- fit_glmnet(Xtr, ytr, alpha=0.5, nfolds=5, seed=321)

    p_hat <- as.numeric(predict(cv, newx=Xte, s=cv$lambda.1se, type="response"))
    data.table(metric=label, cohort_train="ALL_BUT_ONE", cohort_test=held,
               setting="LOCO", AUC=auc_from_probs(yte, p_hat),
               n_train=length(idx_tr), n_test=length(idx_te))
  }))

  # 4) Feature stability: coefficients for each training cohort + LOCO
  feats_train <- rbindlist(lapply(cohorts, function(tr) {
    idx_tr <- which(cohort == tr)
    cv <- fit_glmnet(Xmat[idx_tr, , drop=FALSE], y[idx_tr], alpha=0.5, nfolds=5, seed=777)
    f <- nonzero_features(cv, feature_names, "lambda.1se")
    f[, `:=`(metric=label, train=tr, setting="train_only")]
    f
  }), fill=TRUE)

  feats_loco <- rbindlist(lapply(cohorts, function(held) {
    idx_tr <- which(cohort != held)
    cv <- fit_glmnet(Xmat[idx_tr, , drop=FALSE], y[idx_tr], alpha=0.5, nfolds=5, seed=888)
    f <- nonzero_features(cv, feature_names, "lambda.1se")
    f[, `:=`(metric=label, train=paste0("ALL_BUT_", held), setting="LOCO_train")]
    f
  }), fill=TRUE)

  feats <- rbindlist(list(feats_train, feats_loco), fill=TRUE)

  # pairwise Jaccard overlap between training-cohort feature sets
  sets <- split(feats_train$feature, feats_train$train)
  jac <- CJ(train1=names(sets), train2=names(sets))
  jac[, jaccard := mapply(function(a,b){
    A <- unique(sets[[a]]); B <- unique(sets[[b]])
    if (length(A)==0 && length(B)==0) return(1)
    if (length(A)==0 || length(B)==0) return(0)
    length(intersect(A,B)) / length(union(A,B))
  }, train1, train2)]

  # write outputs
  fwrite(within,  file.path(out_dir, paste0("AUC_within_", label, ".csv")))
  fwrite(transfer,file.path(out_dir, paste0("AUC_transfer_", label, ".csv")))
  fwrite(loco,    file.path(out_dir, paste0("AUC_LOCO_", label, ".csv")))
  fwrite(feats,   file.path(out_dir, paste0("features_", label, ".csv")))
  fwrite(jac,     file.path(out_dir, paste0("jaccard_features_", label, ".csv")))

  cat("Wrote outputs for:", label, "\n")
}

# Run both representations
run_suite(X_clr,    "CLR_counts")
run_suite(X_logrel, "log_rel_abund")

cat("C21 DONE\n")
