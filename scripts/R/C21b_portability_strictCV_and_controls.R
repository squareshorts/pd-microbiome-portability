cat("C21b STARTED\n")

library(data.table)

need_pkgs <- c("glmnet","pROC")
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) {
    stop("Missing package: ", p)
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

Xal <- X[meta[, .(SampleID)], on="SampleID", nomatch=0]
stopifnot(nrow(Xal) == nrow(meta))

y <- as.integer(meta$PD)
if (any(is.na(y))) {
  y <- ifelse(meta$PD %in% c("PD","1",1), 1L, 0L)
}
stopifnot(all(y %in% c(0L,1L)))

cohort <- meta$Cohort

genus_cols <- setdiff(names(Xal), "SampleID")

X_dt <- copy(Xal[, ..genus_cols])
for (col in genus_cols) {
  X_dt[[col]] <- suppressWarnings(as.numeric(X_dt[[col]]))
}
X_dt[is.na(X_dt)] <- 0

X_counts <- as.matrix(X_dt)
storage.mode(X_counts) <- "numeric"

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

X_clr    <- clr_transform(X_counts)
X_logrel <- log(rel_abund(X_counts) + 1e-6)

auc_from_probs <- function(y_true, p_hat) {
  as.numeric(pROC::auc(pROC::roc(y_true, p_hat, quiet=TRUE)))
}

fit_glmnet_cv <- function(Xtr, ytr, alpha=0.5, nfolds=5, seed=1) {
  set.seed(seed)
  cv.glmnet(
    x = Xtr,
    y = ytr,
    family = "binomial",
    alpha = alpha,
    nfolds = nfolds,
    standardize = TRUE
  )
}

# ===========================
# NEW: Stratified fold maker
# ===========================
make_stratified_folds <- function(y, K=5, seed=1) {
  set.seed(seed)
  foldid <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    foldid[idx] <- sample(rep(1:K, length.out=length(idx)))
  }
  foldid
}

# Strict out-of-fold predictions inside cohort (STRATIFIED)
oof_auc <- function(Xc, yc, alpha=0.5, K=5, seed=1) {

  n <- length(yc)
  foldid <- make_stratified_folds(yc, K=K, seed=seed)
  p_hat <- rep(NA_real_, n)

  for (k in 1:K) {
    tr <- which(foldid != k)
    te <- which(foldid == k)

    cv <- fit_glmnet_cv(
      Xc[tr,,drop=FALSE],
      yc[tr],
      alpha=alpha,
      nfolds=K,
      seed=seed+100+k
    )

    p_hat[te] <- as.numeric(
      predict(cv,
              newx=Xc[te,,drop=FALSE],
              s=cv$lambda.1se,
              type="response")
    )
  }

  auc_from_probs(yc, p_hat)
}

nonzero_features <- function(cvfit, feature_names) {
  b <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))
  b <- b[-1, , drop=FALSE]
  nz <- which(abs(b[,1]) > 0)
  data.table(feature=feature_names[nz], coef=b[nz,1], sign=sign(b[nz,1]))
}

out_dir <- "results/paper2_portability_strict"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

run_suite <- function(Xmat, label, perm_reps=50) {

  feature_names <- colnames(Xmat)
  cohorts <- sort(unique(cohort))

  # 1) Strict within-cohort OOF AUC
  within <- rbindlist(lapply(cohorts, function(coh) {
    idx <- which(cohort == coh)
    Xc <- Xmat[idx,,drop=FALSE]
    yc <- y[idx]

    reps <- 20
    aucs <- sapply(1:reps, function(r)
      oof_auc(Xc, yc, alpha=0.5, K=5, seed=1000 + 37*r))

    data.table(metric=label,
               cohort=coh,
               setting="within_strict_oof",
               AUC=mean(aucs),
               AUC_sd=sd(aucs),
               n=length(idx))
  }))

  # 2) Transfer
  transfer <- rbindlist(lapply(cohorts, function(tr) {
    idx_tr <- which(cohort == tr)
    cv <- fit_glmnet_cv(Xmat[idx_tr,,drop=FALSE], y[idx_tr], alpha=0.5, nfolds=5, seed=123)

    rbindlist(lapply(cohorts, function(te) {
      idx_te <- which(cohort == te)
      p_hat <- as.numeric(
        predict(cv,
                newx=Xmat[idx_te,,drop=FALSE],
                s=cv$lambda.1se,
                type="response")
      )
      data.table(metric=label,
                 cohort_train=tr,
                 cohort_test=te,
                 setting="transfer",
                 AUC=auc_from_probs(y[idx_te], p_hat),
                 n_train=length(idx_tr),
                 n_test=length(idx_te))
    }))
  }))

  # 3) LOCO
  loco <- rbindlist(lapply(cohorts, function(held) {
    idx_te <- which(cohort == held)
    idx_tr <- which(cohort != held)

    cv <- fit_glmnet_cv(Xmat[idx_tr,,drop=FALSE], y[idx_tr], alpha=0.5, nfolds=5, seed=321)
    p_hat <- as.numeric(
      predict(cv,
              newx=Xmat[idx_te,,drop=FALSE],
              s=cv$lambda.1se,
              type="response")
    )

    data.table(metric=label,
               cohort_test=held,
               setting="LOCO",
               AUC=auc_from_probs(y[idx_te], p_hat),
               n_train=length(idx_tr),
               n_test=length(idx_te))
  }))

  fwrite(within,   file.path(out_dir, paste0("within_strict_", label, ".csv")))
  fwrite(transfer, file.path(out_dir, paste0("transfer_", label, ".csv")))
  fwrite(loco,     file.path(out_dir, paste0("loco_", label, ".csv")))

  cat("C21b wrote outputs for:", label, "\n")
}

run_suite(X_clr,    "CLR_counts",    perm_reps=50)
run_suite(X_logrel, "log_rel_abund", perm_reps=50)

cat("C21b DONE\n")
