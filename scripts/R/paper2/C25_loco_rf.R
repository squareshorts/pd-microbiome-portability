# scripts/R/paper2/C25_loco_rf.R
cat("C25 LOCO RF STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(ranger)
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
X <- as.data.frame(dt[, ..feat_cols])
y <- as.integer(meta$PD)
stopifnot(all(y %in% c(0,1)))

# Metrics helper (no NA)
calc_metrics <- function(y_true, p_hat, thr=0.5) {
  pred <- ifelse(p_hat >= thr, 1L, 0L)
  tp <- sum(pred==1 & y_true==1)
  tn <- sum(pred==0 & y_true==0)
  fp <- sum(pred==1 & y_true==0)
  fn <- sum(pred==0 & y_true==1)

  sens <- if ((tp+fn) > 0) tp/(tp+fn) else 0
  spec <- if ((tn+fp) > 0) tn/(tn+fp) else 0
  balacc <- (sens + spec)/2

  precision <- if ((tp+fp) > 0) tp/(tp+fp) else 0
  recall <- sens
  f1 <- if ((precision + recall) > 0) 2*precision*recall/(precision+recall) else 0

  list(tp=tp, tn=tn, fp=fp, fn=fn, sens=sens, spec=spec, balacc=balacc, f1=f1)
}

# Output dirs
dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
dir.create("results/paper2/preds", recursive=TRUE, showWarnings=FALSE)

cohorts <- sort(unique(meta$Cohort))
fold_rows <- list()

for (coh in cohorts) {
  cat("\n--- LOCO fold, held-out cohort:", coh, "\n")

  test_idx  <- which(meta$Cohort == coh)
  train_idx <- which(meta$Cohort != coh)

  Xtr <- X[train_idx, , drop=FALSE]
  ytr <- y[train_idx]
  Xte <- X[test_idx, , drop=FALSE]
  yte <- y[test_idx]

  # ranger expects factor outcome for probability
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

  # Predict prob of class "1"
  pred <- predict(fit, data = Xte)$predictions
  p_hat <- as.numeric(pred[, "1"])

  # AUROC
  roc_obj <- tryCatch(pROC::roc(yte, p_hat, quiet=TRUE), error=function(e) NULL)
  auroc <- ifelse(is.null(roc_obj), NA_real_, as.numeric(pROC::auc(roc_obj)))

  # AUPRC (same trapezoid convention used in C24)
  ord <- order(p_hat, decreasing=TRUE)
  yt <- yte[ord]
  tp_cum <- cumsum(yt==1)
  fp_cum <- cumsum(yt==0)
  prec <- tp_cum / pmax(tp_cum + fp_cum, 1)
  rec  <- tp_cum / max(sum(yt==1), 1)
  auprc <- sum(diff(c(0, rec)) * prec)

  m <- calc_metrics(yte, p_hat, thr=0.5)

  out_pred <- data.table(
    SampleID = dt$SampleID[test_idx],
    Cohort   = meta$Cohort[test_idx],
    PD       = yte,
    p_PD     = p_hat
  )
  fwrite(out_pred, file.path("results/paper2/preds", paste0("loco_rf_pred_", coh, ".csv")))

  fold_rows[[length(fold_rows)+1]] <- data.table(
    model="rf_ranger",
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
    tp=m$tp, tn=m$tn, fp=m$fp, fn=m$fn
  )
}

fold_tab <- rbindlist(fold_rows, use.names=TRUE, fill=TRUE)
fwrite(fold_tab, "results/paper2/table_loco_rf_summary.csv")

cat("\nSaved: results/paper2/table_loco_rf_summary.csv\n")
cat("C25 LOCO RF DONE\n")