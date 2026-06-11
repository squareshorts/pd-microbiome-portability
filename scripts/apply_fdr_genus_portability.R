# apply_fdr_genus_portability.R
# Genus-level inferential tests for JMM non-portability manuscript.
# Three hypothesis families, each corrected with BH FDR independently:
#   1. PD main effect      : CLR ~ Cohort + PD
#   2. PD-by-Cohort        : CLR ~ Cohort * PD  (joint F-test vs additive)
#   3. Cross-cohort hetero : per-cohort lm(CLR ~ PD), FE meta + Cochran Q
# Input : results/portability_analysis/dataset_clr.csv  (62 CLR genera)
# Output: results/fdr/*.csv
# Requires: base R only (no external packages)

cat("=== apply_fdr_genus_portability.R ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- Load data ---------------------------------------------------------------
clr_path <- "results/portability_analysis/dataset_clr.csv"
stopifnot(file.exists(clr_path))

dt        <- read.csv(clr_path, stringsAsFactors = FALSE, check.names = FALSE)
feat_cols <- setdiff(names(dt), c("SampleID", "PD", "Cohort"))
n_genera  <- length(feat_cols)

cat("Genera tested:", n_genera, "\n")
cat("Total samples:", nrow(dt), "\n")
cat("Cohort x PD:\n"); print(table(dt$Cohort, dt$PD))

dt$PD     <- factor(dt$PD,     levels = c("0","1"))
dt$Cohort <- factor(dt$Cohort)
dt$Cohort <- relevel(dt$Cohort, ref = "Finland")

X_clr <- as.matrix(dt[, feat_cols])
storage.mode(X_clr) <- "double"

make_row <- function(...) data.frame(..., stringsAsFactors = FALSE)

# ==============================================================================
# FAMILY 1: PD main effect  (CLR ~ Cohort + PD)
# ==============================================================================
cat("\n--- Family 1: PD main effect ---\n")

main_list <- lapply(seq_along(feat_cols), function(j) {
  g   <- feat_cols[j]
  y   <- X_clr[, j]
  fit <- lm(y ~ Cohort + PD, data = dt)
  co  <- summary(fit)$coefficients
  if (!"PD1" %in% rownames(co))
    return(make_row(Genus=g, beta_PD=NA, se_PD=NA, t_PD=NA, p_PD=NA))
  make_row(Genus=g,
           beta_PD=co["PD1",1], se_PD=co["PD1",2],
           t_PD=co["PD1",3],    p_PD=co["PD1",4])
})

main_res           <- do.call(rbind, main_list)
main_res$q_PD      <- p.adjust(main_res$p_PD, method="BH")
main_res$sig_q05   <- main_res$q_PD < 0.05
main_res$direction <- ifelse(main_res$beta_PD > 0, "higher_in_PD", "lower_in_PD")
main_res           <- main_res[order(main_res$p_PD), ]

n1_tested   <- nrow(main_res)
n1_nom      <- sum(main_res$p_PD < 0.05, na.rm=TRUE)
n1_fdr      <- sum(main_res$sig_q05, na.rm=TRUE)
min1_p      <- min(main_res$p_PD, na.rm=TRUE)
min1_q      <- min(main_res$q_PD, na.rm=TRUE)
sig1_genera <- main_res$Genus[!is.na(main_res$sig_q05) & main_res$sig_q05]

cat("  N tested:", n1_tested,"\n")
cat("  Nominal p<0.05:", n1_nom,"\n")
cat("  BH q<0.05:", n1_fdr,"\n")
cat("  Min p:", formatC(min1_p, format="g", digits=4),"\n")
cat("  Min q:", formatC(min1_q, format="g", digits=4),"\n")
cat("\n  Top 10:\n")
print(head(main_res[,c("Genus","beta_PD","se_PD","p_PD","q_PD","direction")],10), row.names=FALSE)

# ==============================================================================
# FAMILY 2: PD-by-Cohort interaction (joint F-test)
# ==============================================================================
cat("\n--- Family 2: PD x Cohort interaction ---\n")

int_list <- lapply(seq_along(feat_cols), function(j) {
  g        <- feat_cols[j]
  y        <- X_clr[, j]
  fit_add  <- lm(y ~ Cohort + PD, data=dt)
  fit_full <- lm(y ~ Cohort * PD, data=dt)
  av       <- anova(fit_add, fit_full)
  p_joint  <- av[["Pr(>F)"]][2]
  co <- summary(fit_full)$coefficients
  get_cf <- function(term) {
    if (term %in% rownames(co)) list(beta=co[term,1], se=co[term,2], p=co[term,4])
    else list(beta=NA, se=NA, p=NA)
  }
  pd_fin <- get_cf("PD1")
  pd_mys <- get_cf("CohortMalaysia:PD1")
  pd_usa <- get_cf("CohortUSA:PD1")
  make_row(Genus=g,
           p_interact_joint=p_joint,
           beta_PD_Finland=pd_fin$beta,
           beta_PD_Malaysia=pd_fin$beta + pd_mys$beta,
           beta_PD_USA=pd_fin$beta + pd_usa$beta,
           beta_PDxMalaysia=pd_mys$beta, se_PDxMalaysia=pd_mys$se, p_PDxMalaysia=pd_mys$p,
           beta_PDxUSA=pd_usa$beta, se_PDxUSA=pd_usa$se, p_PDxUSA=pd_usa$p)
})

int_res                  <- do.call(rbind, int_list)
int_res$q_interact_joint <- p.adjust(int_res$p_interact_joint, method="BH")
int_res$sig_q05          <- int_res$q_interact_joint < 0.05
int_res                  <- int_res[order(int_res$p_interact_joint), ]

n2_tested   <- nrow(int_res)
n2_nom      <- sum(int_res$p_interact_joint < 0.05, na.rm=TRUE)
n2_fdr      <- sum(int_res$sig_q05, na.rm=TRUE)
min2_p      <- min(int_res$p_interact_joint, na.rm=TRUE)
min2_q      <- min(int_res$q_interact_joint, na.rm=TRUE)
sig2_genera <- int_res$Genus[!is.na(int_res$sig_q05) & int_res$sig_q05]

cat("  N tested:", n2_tested,"\n")
cat("  Nominal p<0.05:", n2_nom,"\n")
cat("  BH q<0.05:", n2_fdr,"\n")
cat("  Min p:", formatC(min2_p, format="g", digits=4),"\n")
cat("  Min q:", formatC(min2_q, format="g", digits=4),"\n")
cat("\n  Top 10:\n")
print(head(int_res[,c("Genus","p_interact_joint","q_interact_joint",
                       "beta_PD_Finland","beta_PD_Malaysia","beta_PD_USA")],10), row.names=FALSE)

# ==============================================================================
# FAMILY 3: Cross-cohort heterogeneity (FE meta + Cochran Q)
# ==============================================================================
cat("\n--- Family 3: Cross-cohort heterogeneity ---\n")

cohorts <- sort(unique(as.character(dt$Cohort)))

per_cohort_list <- lapply(cohorts, function(coh) {
  idx  <- which(as.character(dt$Cohort) == coh)
  dt_c <- dt[idx,,drop=FALSE]
  Xc   <- X_clr[idx,,drop=FALSE]
  rows <- lapply(seq_along(feat_cols), function(j) {
    g   <- feat_cols[j]
    y   <- Xc[,j]
    fit <- lm(y ~ PD, data=dt_c)
    co  <- summary(fit)$coefficients
    if (!"PD1" %in% rownames(co)) return(NULL)
    make_row(Cohort=coh, Genus=g,
             beta_PD=co["PD1",1], se_PD=co["PD1",2],
             p_PD=co["PD1",4], n_coh=nrow(dt_c))
  })
  do.call(rbind, rows)
})
per_cohort_res <- do.call(rbind, per_cohort_list)

genus_list <- split(per_cohort_res, per_cohort_res$Genus)
het_list   <- lapply(names(genus_list), function(g) {
  r <- genus_list[[g]]
  r <- r[!is.na(r$se_PD) & r$se_PD > 0, ]
  if (nrow(r) < 2) return(NULL)
  w     <- 1/r$se_PD^2
  bfe   <- sum(w*r$beta_PD)/sum(w)
  se_fe <- sqrt(1/sum(w))
  z     <- bfe/se_fe
  p_fe  <- 2*pnorm(-abs(z))
  Q     <- sum(w*(r$beta_PD - bfe)^2)
  df_q  <- nrow(r)-1
  p_Q   <- pchisq(Q, df=df_q, lower.tail=FALSE)
  I2    <- if (Q > df_q) (Q-df_q)/Q else 0
  make_row(Genus=g, beta_FE=bfe, se_FE=se_fe, z=z, p_FE=p_fe,
           Q=Q, df_Q=df_q, p_Q=p_Q, I2=I2, k=nrow(r))
})
het_res <- do.call(rbind, het_list)

het_res$q_FE      <- p.adjust(het_res$p_FE, method="BH")
het_res$q_Q       <- p.adjust(het_res$p_Q,  method="BH")
het_res$sig_FE    <- het_res$q_FE < 0.05
het_res$sig_Q     <- het_res$q_Q  < 0.05
het_res$direction <- ifelse(het_res$beta_FE > 0, "higher_in_PD", "lower_in_PD")
het_res           <- het_res[order(het_res$p_FE), ]

n3_tested  <- nrow(het_res)
n3_nom_FE  <- sum(het_res$p_FE < 0.05, na.rm=TRUE)
n3_fdr_FE  <- sum(het_res$sig_FE, na.rm=TRUE)
n3_nom_Q   <- sum(het_res$p_Q  < 0.05, na.rm=TRUE)
n3_fdr_Q   <- sum(het_res$sig_Q,  na.rm=TRUE)
min3_p     <- min(het_res$p_FE, na.rm=TRUE)
min3_q     <- min(het_res$q_FE, na.rm=TRUE)
sig3_FE    <- het_res$Genus[!is.na(het_res$sig_FE) & het_res$sig_FE]
sig3_Q     <- het_res$Genus[!is.na(het_res$sig_Q)  & het_res$sig_Q]

cat("  N tested:", n3_tested,"\n")
cat("  Nominal p_FE<0.05:", n3_nom_FE,"\n")
cat("  BH q_FE<0.05:", n3_fdr_FE,"\n")
cat("  Nominal p_Q<0.05:", n3_nom_Q,"\n")
cat("  BH q_Q<0.05:", n3_fdr_Q,"\n")
cat("\n  Top 10:\n")
print(head(het_res[,c("Genus","beta_FE","se_FE","p_FE","q_FE","Q","p_Q","I2","direction")],10), row.names=FALSE)

# ==============================================================================
# PC1 alignment (cohort-level, BH pre-applied in C33)
# ==============================================================================
pc_align_path <- "results/paper2/table_pd_effect_alignment_angles.csv"
pc_align_out  <- NULL
if (file.exists(pc_align_path)) {
  pc <- read.csv(pc_align_path, stringsAsFactors=FALSE)
  pc_align_out <- data.frame(
    cohort=pc$cohort, n=pc$n, prev=pc$prev,
    abs_cos_pc1=pc$abs_cos_pc1, angle_pc1_deg=pc$angle_pc1_deg,
    p_perm_align_pc1=pc$p_perm_align_pc1,
    q_perm_align_pc1=pc$p_perm_align_pc1_adj,
    sig_q05=pc$p_perm_align_pc1_adj < 0.05,
    stringsAsFactors=FALSE)
  cat("\n--- PC1 alignment (cohort-level) ---\n")
  print(pc_align_out, row.names=FALSE)
}

# ==============================================================================
# Save outputs
# ==============================================================================
dir.create("results/fdr", showWarnings=FALSE, recursive=TRUE)

write.csv(main_res,       "results/fdr/genus_pd_main_effects_fdr.csv",                row.names=FALSE)
write.csv(int_res,        "results/fdr/genus_disease_by_cohort_interactions_fdr.csv", row.names=FALSE)
write.csv(per_cohort_res, "results/fdr/genus_per_cohort_pd_effects.csv",              row.names=FALSE)
write.csv(het_res,        "results/fdr/genus_heterogeneity_fdr.csv",                  row.names=FALSE)
if (!is.null(pc_align_out))
  write.csv(pc_align_out, "results/fdr/cohort_pc_alignment_fdr.csv",                  row.names=FALSE)

# ==============================================================================
# FDR summary
# ==============================================================================
now_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
fdr_rows <- list(
  make_row(family="PD_main_effect", model="CLR ~ Cohort + PD",
           input_file=clr_path, p_column="p_PD", q_column="q_PD",
           n_tested=n1_tested, n_nominal_p05=n1_nom, n_fdr_q05=n1_fdr,
           min_p=min1_p, min_q=min1_q,
           sig_genera=paste(sig1_genera, collapse="; "),
           generated_at=now_str, script="scripts/apply_fdr_genus_portability.R"),
  make_row(family="PD_by_Cohort_interaction",
           model="CLR ~ Cohort * PD (joint F-test vs additive)",
           input_file=clr_path, p_column="p_interact_joint", q_column="q_interact_joint",
           n_tested=n2_tested, n_nominal_p05=n2_nom, n_fdr_q05=n2_fdr,
           min_p=min2_p, min_q=min2_q,
           sig_genera=paste(sig2_genera, collapse="; "),
           generated_at=now_str, script="scripts/apply_fdr_genus_portability.R"),
  make_row(family="Cross_cohort_heterogeneity_FE",
           model="Per-cohort lm(CLR ~ PD); FE meta + Cochran Q",
           input_file=clr_path, p_column="p_FE", q_column="q_FE",
           n_tested=n3_tested, n_nominal_p05=n3_nom_FE, n_fdr_q05=n3_fdr_FE,
           min_p=min3_p, min_q=min3_q,
           sig_genera=paste(sig3_FE, collapse="; "),
           generated_at=now_str, script="scripts/apply_fdr_genus_portability.R"))

if (!is.null(pc_align_out)) {
  fdr_rows <- c(fdr_rows, list(
    make_row(family="PC1_alignment_permutation",
             model="Permutation B=2000: disease-vector vs PC1, per cohort",
             input_file=pc_align_path,
             p_column="p_perm_align_pc1",
             q_column="q_perm_align_pc1 (BH; source C33)",
             n_tested=nrow(pc_align_out),
             n_nominal_p05=sum(pc_align_out$p_perm_align_pc1<0.05, na.rm=TRUE),
             n_fdr_q05=sum(pc_align_out$sig_q05, na.rm=TRUE),
             min_p=min(pc_align_out$p_perm_align_pc1, na.rm=TRUE),
             min_q=min(pc_align_out$q_perm_align_pc1, na.rm=TRUE),
             sig_genera="(cohort-level; not genus-level)",
             generated_at=now_str,
             script="scripts/R/paper2/C33_pd_effect_alignment_angles.R")))
}
fdr_summary <- do.call(rbind, fdr_rows)
write.csv(fdr_summary, "results/fdr/fdr_summary.csv", row.names=FALSE)

# ==============================================================================
# Console report
# ==============================================================================
cat("\n\n=== FINAL FDR SUMMARY ===\n")
for (i in seq_len(nrow(fdr_summary))) {
  r <- fdr_summary[i,]
  cat(sprintf("\n[%s]\n", r$family))
  cat(sprintf("  Model:          %s\n", r$model))
  cat(sprintf("  N tested:       %d\n", r$n_tested))
  cat(sprintf("  Nominal p<0.05: %d\n", r$n_nominal_p05))
  cat(sprintf("  BH q<0.05:      %d\n", r$n_fdr_q05))
  cat(sprintf("  Min p: %.4g | Min q: %.4g\n", r$min_p, r$min_q))
  if (r$n_fdr_q05 > 0 && !grepl("cohort-level", r$sig_genera))
    cat(sprintf("  Significant genera: %s\n", r$sig_genera))
  else if (r$n_fdr_q05 == 0)
    cat("  No genera reached q<0.05\n")
}
cat("\n=== HETEROGENEITY DETAIL ===\n")
cat("Genera sig heterogeneity (q_Q<0.05):", n3_fdr_Q, "\n")
if (n3_fdr_Q > 0) cat("Names:", paste(sig3_Q, collapse="; "), "\n")
cat("\n=== TOP RESULTS ===\n")
cat("\nPD main effect top 3:\n")
print(head(main_res[,c("Genus","beta_PD","se_PD","p_PD","q_PD","direction")],3), row.names=FALSE)
cat("\nInteraction top 3:\n")
print(head(int_res[,c("Genus","p_interact_joint","q_interact_joint",
                       "beta_PD_Finland","beta_PD_Malaysia","beta_PD_USA")],3), row.names=FALSE)
cat("\nHeterogeneity FE top 3:\n")
print(head(het_res[,c("Genus","beta_FE","se_FE","p_FE","q_FE","I2","direction")],3), row.names=FALSE)
cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
