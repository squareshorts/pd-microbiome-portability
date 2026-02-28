library(data.table)
library(ggplot2)

base <- "results/paper2_portability_strict"
outdir <- "figures"
dir.create(outdir, showWarnings = FALSE)

# =========================
# FIGURE 1 — Performance
# =========================

within_clr <- fread(file.path(base,"within_strict_CLR_counts.csv"))
loco_clr   <- fread(file.path(base,"loco_CLR_counts.csv"))

within_rel <- fread(file.path(base,"within_strict_log_rel_abund.csv"))
loco_rel   <- fread(file.path(base,"loco_log_rel_abund.csv"))

perf_clr <- rbindlist(list(
  within_clr[,.(cohort, AUC, type="Within")],
  loco_clr[,.(cohort=cohort_test, AUC, type="LOCO")]
))

perf_rel <- rbindlist(list(
  within_rel[,.(cohort, AUC, type="Within")],
  loco_rel[,.(cohort=cohort_test, AUC, type="LOCO")]
))

perf_clr[, cohort := factor(cohort, levels=c("finland","malaysia","usa"))]
perf_rel[, cohort := factor(cohort, levels=c("finland","malaysia","usa"))]

p1 <- ggplot(perf_clr, aes(x=cohort, y=AUC, fill=type)) +
  geom_col(position=position_dodge(width=0.8), width=0.7) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  theme_bw()

p2 <- ggplot(perf_rel, aes(x=cohort, y=AUC, fill=type)) +
  geom_col(position=position_dodge(width=0.8), width=0.7) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  theme_bw()

ggsave(file.path(outdir,"Figure1_CLR_performance.pdf"), p1, width=6, height=4)
ggsave(file.path(outdir,"Figure1_logrel_performance.pdf"), p2, width=6, height=4)

# =========================
# FIGURE 2 — Transfer Matrix
# =========================

transfer_clr <- fread(file.path(base,"transfer_CLR_counts.csv"))

mat_clr <- dcast(
  transfer_clr,
  cohort_train ~ cohort_test,
  value.var="AUC"
)

m_clr <- melt(mat_clr,
              id.vars="cohort_train",
              variable.name="cohort_test",
              value.name="AUC")

p_heat <- ggplot(m_clr,
                 aes(x=cohort_train,
                     y=cohort_test,
                     fill=AUC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(0.45,1)) +
  theme_bw()

ggsave(file.path(outdir,"Figure2_transfer_CLR_heatmap.pdf"),
       p_heat, width=5, height=5)

# =========================
# FIGURE 3 — Jaccard
# =========================

jac <- fread(file.path(base,"jaccard_trainonly_CLR_counts.csv"))

p_jac <- ggplot(jac,
                aes(x=train1,
                    y=train2,
                    fill=jaccard)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(0,1)) +
  theme_bw()

ggsave(file.path(outdir,"Figure3_jaccard_CLR.pdf"),
       p_jac, width=5, height=5)

# =========================
# FIGURE 4 — Top Coefficients
# =========================

feat <- fread(file.path(base,"features_trainonly_CLR_counts.csv"))
feat[,abscoef:=abs(coef)]

top_feats <- feat[
  order(-abscoef),
  head(.SD,10),
  by=train
]

p_coef <- ggplot(top_feats,
                 aes(x=reorder(feature,abscoef),
                     y=coef,
                     fill=train)) +
  geom_col() +
  coord_flip() +
  theme_bw()

ggsave(file.path(outdir,"Figure4_top_features_CLR.pdf"),
       p_coef, width=7, height=6)

cat("Figures rebuilt successfully.\n")
