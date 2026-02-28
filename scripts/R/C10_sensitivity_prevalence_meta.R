library(data.table)

meta_file  <- "data/processed/ipd_sex/ipd_metadata.csv"
count_file <- "data/processed/ipd_sex/ipd_genus_counts_union.csv"

if (!file.exists(meta_file))  stop("Missing metadata file.")
if (!file.exists(count_file)) stop("Missing genus count file.")

meta <- fread(meta_file)
X    <- fread(count_file)

setkey(meta, SampleID)
setkey(X, SampleID)

dt <- X[meta]

if (nrow(dt) == 0)
  stop("Join failed: no overlapping SampleID.")

required_cols <- c("SampleID","PD","Sex","Cohort")
if (!all(required_cols %in% names(dt)))
  stop("Missing required columns in merged data.")

genus_cols <- setdiff(names(dt), required_cols)

clr_transform <- function(mat) {
  mat_pc <- mat + 1
  log_mat <- log(mat_pc)
  log_mat - rowMeans(log_mat)
}

run_meta <- function(prev_thr) {

  cat("\nRunning prevalence threshold:", prev_thr, "\n")

  prev <- colMeans(as.matrix(dt[, ..genus_cols]) > 0)
  keep <- genus_cols[prev >= prev_thr]

  cat("Genera retained:", length(keep), "\n")

  if (length(keep) == 0)
    return(NULL)

  cohorts <- sort(unique(dt$Cohort))

  cohort_res <- rbindlist(lapply(cohorts, function(coh) {

    d <- dt[Cohort == coh]

    if (nrow(d) < 20)
      return(NULL)

    mat <- as.matrix(d[, ..keep])
    storage.mode(mat) <- "numeric"

    clr <- clr_transform(mat)

    df <- data.table(
      PD  = factor(d$PD),
      Sex = factor(d$Sex)
    )

    if (!("0" %in% levels(df$PD))) return(NULL)
    if (!("Female" %in% levels(df$Sex))) return(NULL)

    df[, PD  := relevel(PD, ref = "0")]
    df[, Sex := relevel(Sex, ref = "Female")]

    rbindlist(lapply(seq_along(keep), function(j) {

      g <- keep[j]
      y <- clr[, j]

      fit <- try(lm(y ~ PD * Sex, data = df), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)

      coefs <- summary(fit)$coefficients

      if (!("PD1:SexMale" %in% rownames(coefs)))
        return(NULL)

      data.table(
        Cohort = coh,
        Genus  = g,
        beta   = coefs["PD1:SexMale",1],
        se     = coefs["PD1:SexMale",2]
      )
    }))

  }), fill = TRUE)

  if (is.null(cohort_res) || nrow(cohort_res) == 0)
    return(NULL)

  meta_res <- cohort_res[, {

    if (.N < 2)
      return(NULL)

    w <- 1 / (se^2)

    beta_fe <- sum(w * beta) / sum(w)
    se_fe   <- sqrt(1 / sum(w))

    z    <- beta_fe / se_fe
    p_fe <- 2 * pnorm(-abs(z))

    Q  <- sum(w * (beta - beta_fe)^2)
    df <- .N - 1
    p_Q <- pchisq(Q, df = df, lower.tail = FALSE)

    I2 <- ifelse(Q > df, (Q - df) / Q, 0)

    .(
      beta_FE = beta_fe,
      se_FE   = se_fe,
      z       = z,
      p_FE    = p_fe,
      Q       = Q,
      df      = df,
      p_Q     = p_Q,
      I2      = I2,
      k       = .N
    )

  }, by = Genus]

  if (is.null(meta_res) || nrow(meta_res) == 0)
    return(NULL)

  meta_res[, q_FE := p.adjust(p_FE, method = "BH")]
  setorder(meta_res, p_FE)

  meta_res
}

dir.create("results/sex_paper/sensitivity", showWarnings = FALSE, recursive = TRUE)

for (thr in c(0.05, 0.10, 0.20)) {

  res <- run_meta(thr)

  if (is.null(res)) {
    cat("No valid genera retained at threshold", thr, "\n")
    next
  }

  out_file <- sprintf(
    "results/sex_paper/sensitivity/meta_prev_%02d.csv",
    as.integer(thr * 100)
  )

  fwrite(res, out_file)

  cat("Wrote:", out_file, "\n")
  print(res[1:min(10,.N),
            .(Genus, p_FE, q_FE, beta_FE, se_FE, I2, k)])
}
