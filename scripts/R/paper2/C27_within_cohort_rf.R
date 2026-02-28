# scripts/R/paper2/C27_within_cohort_rf.R
cat("C27 WITHIN-COHORT RF STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(pROC)
})

x_path    <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"
stopifnot(file.exists(x_path), file.exists(meta_path))

Xdt  <- fread(x_path); setnames(Xdt, 1, "SampleID")
meta <- fread(meta_path)
stopifnot(all(c("SampleID","PD","Cohort") %in% names(meta)))

setkey(meta, SampleID)
setkey(Xdt, SampleID)
dt <- Xdt[meta[,.(SampleID)], on="SampleID", nomatch=0]
stopifnot(nrow(dt)==nrow(meta))

feat_cols <- setdiff(names(dt), "SampleID")
X_all <- as.data.frame(dt[, ..feat_cols])
y_all <- as.integer(meta$PD)
stopifnot(all(y_all %in% c(0,1)))

calc_metrics <- function(y_true, p_hat, thr=0.5) {
  pred <- ifelse(p_hat >= thr, 1L, 0L)
  tp <- sum(pred==1 & y_true==1)
  tn <- sum(pred==0 & y_true==0)
  fp <- sum(pred==1 & y_true==0)
  fn <- sum(pred==0 & y_true==1)
  sens <- if ((tp+fn) > 0) tp/(tp+fn) else 0
  spec <- if ((tn+fp) > 0) tn/(tn+fp) else 0
  balacc <- (sens+spec)/2
  precision <- if ((tp+fp) > 0) tp/(tp+fp) else 0
  recall <- sens
  f1 <- if ((precision+recall) > 0) 2*precision*recall/(precision+recall) else 0
  list(sens=sens, spec=spec, balacc=balacc, f1=f1)
}

make_strat_folds <- function(y, K=5, seed=1) {
  set.seed(seed)
  idx1 <- which(y==1); idx0 <- which(y==0)
  idx1 <- sample(idx1); idx0 <- sample(idx0)
  folds <- vector("list", K)
  for (k in 1:K) folds[[k]] <- integer(0)
  for (i in seq_along(idx1)) folds[[ (i-1) %% K + 1 ]] <- c(folds[[ (i-1) %% K + 1 ]], idx1[i])
  for (i in seq_along(idx0)) folds[[ (i-1) %% K + 1 ]] <- c(folds[[ (i-1) %% K + 1 ]], idx0[i])
  folds
}

dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
dir.create("results/paper2/within", recursive=TRUE, showWarnings=FALSE)

cohorts <- sort(unique(meta$Cohort))
out_rows <- list()

K <- 5
R <- 10

for (coh in cohorts) {
  cat("\n--- Cohort:", coh, "\n")
  idx <- which(meta$Cohort == coh)
  Xc <- X_all[idx, , drop=FALSE]
  yc <- y_all[idx]
  n <- length(yc)
  prev <- mean(yc)

  fold_metrics <- list()

  for (r in 1:R) {
    folds <- make_strat_folds(yc, K=K, seed=200 + r)
    for (k in 1:K) {
      test_idx <- folds[[k]]
      train_idx <- setdiff(seq_len(n), test_idx)

      Xtr <- Xc[train_idx, , drop=FALSE]
      ytr <- yc[train_idx]
      Xte <- Xc[test_idx, , drop=FALSE]
      yte <- yc[test_idx]

      dtr <- data.table(PD=factor(ytr, levels=c(0,1)), Xtr)

      set.seed(1)
      fit <- ranger(
        formula = PD ~ .,
        data = as.data.frame(dtr),
        probability = TRUE,
        num.trees = 1000,
        mtry = max(1, floor(sqrt(ncol(Xtr)))),
        min.node.size = 5,
        respect.unordered.factors = "order",
        seed = 1
      )

      pred <- predict(fit, data=Xte)$predictions
      p_hat <- as.numeric(pred[, "1"])

      roc_obj <- tryCatch(pROC::roc(yte, p_hat, quiet=TRUE), error=function(e) NULL)
      auroc <- ifelse(is.null(roc_obj), NA_real_, as.numeric(pROC::auc(roc_obj)))

      ord <- order(p_hat, decreasing=TRUE)
      yt <- yte[ord]
      tp_cum <- cumsum(yt==1)
      fp_cum <- cumsum(yt==0)
      prec <- tp_cum / pmax(tp_cum + fp_cum, 1)
      rec  <- tp_cum / max(sum(yt==1), 1)
      auprc <- sum(diff(c(0, rec)) * prec)

      m <- calc_metrics(yte, p_hat, thr=0.5)

      fold_metrics[[length(fold_metrics)+1]] <- data.table(
        cohort=coh, rep_id=r, fold=k,
        n_train=length(train_idx), n_test=length(test_idx),
        prev_test=mean(yte),
        auroc=auroc, auprc=auprc,
        f1=m$f1, balacc=m$balacc, sensitivity=m$sens, specificity=m$spec
      )
    }
  }

  fold_tab <- rbindlist(fold_metrics)
  fwrite(fold_tab, file.path("results/paper2/within", paste0("within_rf_", coh, "_folds.csv")))

  out_rows[[length(out_rows)+1]] <- fold_tab[, .(
    model="within_rf_ranger",
    cohort=coh,
    n=n,
    prev=prev,
    auroc_mean=mean(auroc, na.rm=TRUE),
    auroc_sd=sd(auroc, na.rm=TRUE),
    auprc_mean=mean(auprc, na.rm=TRUE),
    auprc_sd=sd(auprc, na.rm=TRUE),
    f1_mean=mean(f1, na.rm=TRUE),
    f1_sd=sd(f1, na.rm=TRUE),
    balacc_mean=mean(balacc, na.rm=TRUE),
    balacc_sd=sd(balacc, na.rm=TRUE)
  )]
}

sum_tab <- rbindlist(out_rows, use.names=TRUE, fill=TRUE)
fwrite(sum_tab, "results/paper2/table_within_cohort_rf_summary.csv")

cat("\nSaved: results/paper2/table_within_cohort_rf_summary.csv\n")
cat("C27 WITHIN-COHORT RF DONE\n")