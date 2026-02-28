# scripts/R/paper2/C24_loco_logreg.R
cat("C24 LOCO LOGREG STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(glmnet)
})

# Paths (Paper2 only)
x_path    <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"

stopifnot(file.exists(x_path), file.exists(meta_path))

Xdt  <- fread(x_path)
meta <- fread(meta_path)

# Standardize ID column name
setnames(Xdt, 1, "SampleID")
stopifnot("SampleID" %in% names(meta), "PD" %in% names(meta), "Cohort" %in% names(meta))

# Align rows strictly by SampleID
setkey(meta, SampleID)
setkey(Xdt, SampleID)
dt <- Xdt[meta[, .(SampleID)], on="SampleID", nomatch=0]
stopifnot(nrow(dt) == nrow(meta))

# Extract matrices
feat_cols <- setdiff(names(dt), "SampleID")
X <- as.matrix(dt[, ..feat_cols])
y <- as.integer(meta$PD)
stopifnot(all(y %in% c(0,1)))

cohorts <- sort(unique(meta$Cohort))

# Metrics helper
calc_metrics <- function(y_true, p_hat, thr=0.5) {
  pred <- ifelse(p_hat >= thr, 1L, 0L)

  tp <- sum(pred==1 & y_true==1)
  tn <- sum(pred==0 & y_true==0)
  fp <- sum(pred==1 & y_true==0)
  fn <- sum(pred==0 & y_true==1)

  sens <- if ((tp+fn) > 0) tp/(tp+fn) else 0
  spec <- if ((tn+fp) > 0) tn/(tn+fp) else 0
  balacc <- (sens + spec) / 2

  precision <- if ((tp+fp) > 0) tp/(tp+fp) else 0
  recall    <- sens

  f1 <- if ((precision + recall) > 0)
    2 * precision * recall / (precision + recall)
  else
    0

  list(tp=tp, tn=tn, fp=fp, fn=fn,
       sens=sens, spec=spec,
       balacc=balacc, f1=f1)
}

# Output dirs
dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
dir.create("results/paper2/preds", recursive=TRUE, showWarnings=FALSE)

fold_rows <- list()

for (coh in cohorts) {
  cat("\n--- LOCO fold, held-out cohort:", coh, "\n")

  test_idx  <- which(meta$Cohort == coh)
  train_idx <- which(meta$Cohort != coh)

  Xtr <- X[train_idx, , drop=FALSE]
  ytr <- y[train_idx]
  Xte <- X[test_idx, , drop=FALSE]
  yte <- y[test_idx]

  # Regularized logistic regression with internal CV for lambda only (within training data)
  # NOTE: standardize=FALSE because CLR already scaled-ish; keep deterministic.
  set.seed(1)
  cvfit <- cv.glmnet(
    x=Xtr, y=ytr, family="binomial",
    alpha=0,                 # ridge (stable under collinearity)
    nfolds=5,
    type.measure="deviance",
    standardize=FALSE
  )
  lam <- cvfit$lambda.min

  fit <- glmnet(
    x=Xtr, y=ytr, family="binomial",
    alpha=0,
    lambda=lam,
    standardize=FALSE
  )

  # Predict probabilities on held-out cohort
  p_hat <- as.numeric(predict(fit, newx=Xte, type="response"))

  # AUROC
  roc_obj <- tryCatch(pROC::roc(yte, p_hat, quiet=TRUE), error=function(e) NULL)
  auroc <- ifelse(is.null(roc_obj), NA_real_, as.numeric(pROC::auc(roc_obj)))

  # AUPRC (manual via trapezoid on precision-recall)
  # Build PR curve
  ord <- order(p_hat, decreasing=TRUE)
  yt <- yte[ord]
  tp_cum <- cumsum(yt==1)
  fp_cum <- cumsum(yt==0)
  prec <- tp_cum / pmax(tp_cum + fp_cum, 1)
  rec  <- tp_cum / max(sum(yt==1), 1)
  # PR AUC via trapezoid on recall
  auprc <- sum(diff(c(0, rec)) * prec)

  # Thresholded metrics at 0.5
  m <- calc_metrics(yte, p_hat, thr=0.5)

  # Save predictions
  out_pred <- data.table(
    SampleID = dt$SampleID[test_idx],
    Cohort   = meta$Cohort[test_idx],
    PD       = yte,
    p_PD     = p_hat
  )
  fwrite(out_pred, file.path("results/paper2/preds", paste0("loco_logreg_pred_", coh, ".csv")))

  fold_rows[[length(fold_rows)+1]] <- data.table(
    model="logreg_ridge",
    heldout_cohort=coh,
    n_train=length(train_idx),
    n_test=length(test_idx),
    prev_test=mean(yte),
    auroc=auroc,
    auprc=auprc,
    f1=m$f1,
    balacc=m$balacc,
    sensitivity=m$sens,
    specificity=m$spec,
    tp=m$tp, tn=m$tn, fp=m$fp, fn=m$fn,
    lambda=lam
  )
}

fold_tab <- rbindlist(fold_rows, use.names=TRUE, fill=TRUE)
fwrite(fold_tab, "results/paper2/table_loco_logreg_summary.csv")

cat("\nSaved: results/paper2/table_loco_logreg_summary.csv\n")
cat("C24 LOCO LOGREG DONE\n")