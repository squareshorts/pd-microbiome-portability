library(data.table)
library(ggplot2)

in_path <- "results/sex_paper/forest/forest_top10.csv"
if (!file.exists(in_path)) stop("Missing: ", in_path)

d <- fread(in_path)

required_cols <- c("Genus","Cohort","beta","se","beta_FE","se_FE","p_FE","q_FE","I2")
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse=", "))

dir.create("results/sex_paper/forest", showWarnings=FALSE, recursive=TRUE)
pdf("results/sex_paper/forest/forest_top10.pdf", width=7, height=4)

for (g in unique(d$Genus)) {

  dg <- d[Genus == g]

  cohort_rows <- dg[, .(
    Study = Cohort,
    beta  = beta,
    se    = se
  )]

  meta_row <- dg[1, .(
    Study = "Meta (FE)",
    beta  = beta_FE,
    se    = se_FE
  )]

  dd <- rbindlist(list(cohort_rows, meta_row), use.names=TRUE)
  dd <- dd[!is.na(se) & se > 0]

  dd[, lo := beta - 1.96 * se]
  dd[, hi := beta + 1.96 * se]

  study_levels <- c(sort(unique(cohort_rows$Study)), "Meta (FE)")
  dd[, Study := factor(Study, levels = rev(study_levels))]

  subtitle <- dg[1, sprintf("p=%.3g | q=%.3g | I²=%.2f%%",
                            p_FE, q_FE, I2 * 100)]

  p <- ggplot(dd, aes(x = beta, y = Study)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_point(size = 2) +
    labs(title = g,
         subtitle = subtitle,
         x = "PD × Sex interaction (CLR units)",
         y = NULL) +
    theme_bw()

  print(p)
}

dev.off()
cat("Wrote results/sex_paper/forest/forest_top10.pdf\n")
