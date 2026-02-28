library(data.table)

meta <- fread("data/processed/ipd_sex/ipd_metadata.csv")
X <- fread("data/processed/ipd_sex/ipd_genus_counts_union.csv")

setkey(meta, SampleID); setkey(X, SampleID)
dt <- X[meta]

genus_cols <- setdiff(names(dt), c("SampleID","PD","Sex","Cohort"))

# CLR
mat <- as.matrix(dt[, ..genus_cols])
storage.mode(mat) <- "numeric"
mat_pc <- mat + 1
log_mat <- log(mat_pc)
clr_mat <- log_mat - rowMeans(log_mat)

# Design
df <- data.table(
  PD = factor(dt$PD),
  Sex = factor(dt$Sex),
  Cohort = factor(dt$Cohort)
)

# Ensure reference levels are stable
df[, Sex := relevel(Sex, ref="Female")]
df[, PD := relevel(PD, ref="0")]
df[, Cohort := relevel(Cohort, ref="finland")]

results <- rbindlist(lapply(seq_along(genus_cols), function(j){
  g <- genus_cols[j]
  y <- clr_mat[, j]

  fit <- lm(y ~ Cohort + PD*Sex, data=df)
  co <- summary(fit)$coefficients

  # term names:
  # PD1 = PD effect in females
  # SexMale = sex main effect among controls
  # PD1:SexMale = interaction (extra PD effect in males vs females)

  get_term <- function(term) {
    if (term %in% rownames(co)) {
      return(list(beta=co[term,1], se=co[term,2], t=co[term,3], p=co[term,4]))
    } else {
      return(list(beta=NA_real_, se=NA_real_, t=NA_real_, p=NA_real_))
    }
  }

  pd_f <- get_term("PD1")
  int  <- get_term("PD1:SexMale")

  # PD effect in males = PD1 + interaction
  beta_pd_male <- pd_f$beta + int$beta

  data.table(
    Genus = g,
    beta_PD_female = pd_f$beta,
    p_PD_female = pd_f$p,
    beta_interaction = int$beta,
    p_interaction = int$p,
    beta_PD_male = beta_pd_male
  )
}))

results[, q_interaction := p.adjust(p_interaction, method="BH")]
results[, q_PD_female := p.adjust(p_PD_female, method="BH")]

dir.create("results/sex_paper", showWarnings=FALSE, recursive=TRUE)
fwrite(results, "results/sex_paper/genus_interaction_results.csv")

# Print top interaction hits (unadjusted and adjusted)
cat("Top 15 genera by smallest p_interaction:\n")
print(results[order(p_interaction)][1:15,
      .(Genus, beta_interaction, p_interaction, q_interaction, beta_PD_female, beta_PD_male)])
