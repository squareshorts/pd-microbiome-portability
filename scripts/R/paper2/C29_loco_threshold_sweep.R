cat("C29 LOCO THRESHOLD SWEEP STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

x_path    <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"
stopifnot(file.exists(x_path), file.exists(meta_path))

Xdt  <- fread(x_path)
meta <- fread(meta_path)
setnames(Xdt, 1, "SampleID")

setkey(meta, SampleID)
setkey(Xdt, SampleID)
dt <- Xdt[meta[, .(SampleID)], on="SampleID", nomatch=0]
stopifnot(nrow(dt) == nrow(meta))

feat_cols <- setdiff(names(dt), "SampleID")
X <- as.matrix(dt[, ..feat_cols])
y <- as.integer(meta$PD)
cohorts <- sort(unique(meta$Cohort))

calc_basic <- function(y_true, p_hat, thr) {
  pred <- ifelse(p_hat >= thr, 1L, 0L)
  tp <- sum(pred==1 & y_true==1); tn <- sum(pred==0 & y_true==0)
  fp <- sum(pred==1 & y_true==0); fn <- sum(pred==0 & y_true==1)
  sens <- ifelse((tp+fn)>0, tp/(tp+fn), NA_real_)
  spec <- ifelse((tn+fp)>0, tn/(tn+fp), NA_real_)
  prec <- ifelse((tp+fp)>0, tp/(tp+fp), NA_real_)
  f1 <- ifelse(is.na(prec) | is.na(sens) | (prec+sens)==0, NA_real_, 2*prec*sens/(prec+sens))
  balacc <- mean(c(sens,spec), na.rm=TRUE)
  list(tp=tp, tn=tn, fp=fp, fn=fn, sens=sens, spec=spec, prec=prec, f1=f1, balacc=balacc)
}

dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
out_rows <- list()

for (coh in cohorts) {
  cat("\n--- Held-out cohort:", coh, "\n")
  test_idx  <- which(meta$Cohort == coh)
  train_idx <- which(meta$Cohort != coh)

  Xtr <- X[train_idx,,drop=FALSE]; ytr <- y[train_idx]
  Xte <- X[test_idx,,drop=FALSE];  yte <- y[test_idx]

  set.seed(1)
  cvfit <- cv.glmnet(x=Xtr, y=ytr, family="binomial", alpha=0, nfolds=5,
                     type.measure="deviance", standardize=FALSE)
  lam <- cvfit$lambda.min
  fit <- glmnet(x=Xtr, y=ytr, family="binomial", alpha=0, lambda=lam, standardize=FALSE)

  p_tr <- as.numeric(predict(fit, newx=Xtr, type="response"))
  p_te <- as.numeric(predict(fit, newx=Xte, type="response"))

  # Sweep thresholds on TRAIN, pick best by F1 and by Youden J (sens+spec-1)
  thrs <- seq(0.05, 0.95, by=0.01)
  sweep <- rbindlist(lapply(thrs, function(t) {
    m <- calc_basic(ytr, p_tr, t)
    data.table(thr=t, f1=m$f1, youden=(m$sens + m$spec - 1), balacc=m$balacc)
  }))

  best_f1 <- sweep[which.max(f1)][1]
  best_j  <- sweep[which.max(youden)][1]

  # Apply chosen thresholds to TEST
  m05 <- calc_basic(yte, p_te, 0.5)
  mf1 <- calc_basic(yte, p_te, best_f1$thr)
  mj  <- calc_basic(yte, p_te, best_j$thr)

  out_rows[[length(out_rows)+1]] <- data.table(
    cohort=coh,
    thr_train_best_f1 = best_f1$thr,
    thr_train_best_j  = best_j$thr,
    test_f1_thr05     = m05$f1,
    test_f1_best_f1   = mf1$f1,
    test_f1_best_j    = mj$f1,
    test_balacc_thr05   = m05$balacc,
    test_balacc_best_f1 = mf1$balacc,
    test_balacc_best_j  = mj$balacc
  )
}

out <- rbindlist(out_rows, use.names=TRUE, fill=TRUE)
out_path <- "results/paper2/table_loco_threshold_sweep_logreg.csv"
fwrite(out, out_path)
cat("\nSaved:", out_path, "\n")
cat("C29 LOCO THRESHOLD SWEEP DONE\n")