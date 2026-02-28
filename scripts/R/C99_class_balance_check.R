# ============================================
# C99_class_balance_check.R
# Purpose: Inspect PD/control balance per cohort
# ============================================

cat("Script started\n")

meta_path <- "data/processed/ipd_sex/ipd_metadata.csv"
cat("Reading:", meta_path, "\n")

meta <- read.csv(meta_path, stringsAsFactors = FALSE)

cat("Loaded metadata with dimensions:\n")
print(dim(meta))

cat("\nColumn names:\n")
print(colnames(meta))

# Correct column names from your file
cohort_col <- "Cohort"
pd_col     <- "PD"

cat("\nRaw counts (Cohort x PD):\n")
tab <- table(meta[[cohort_col]], meta[[pd_col]])
print(tab)

cat("\nRow proportions (within cohort):\n")
prop_tab <- prop.table(tab, margin = 1)
print(round(prop_tab, 3))

cat("\nTotal subjects per cohort:\n")
print(rowSums(tab))

cat("\nOverall totals:\n")
print(colSums(tab))
