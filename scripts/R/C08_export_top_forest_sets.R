library(data.table)

cat("getwd():", getwd(), "\n")

forest_path <- "results/sex_paper/forest_table_interactions.csv"
meta_path   <- "results/sex_paper/meta_interaction_FE_heterogeneity.csv"

if (!file.exists(forest_path)) stop("Missing: ", forest_path)
if (!file.exists(meta_path))   stop("Missing: ", meta_path)

forest <- fread(forest_path)
meta   <- fread(meta_path)

# Keep only genera present in ≥2 cohorts
meta <- meta[k >= 2]

if (nrow(meta) == 0) stop("No genera with k >= 2.")

# Order by p-value
setorder(meta, p_FE)

topN <- min(10L, nrow(meta))
meta_top <- meta[1:topN, .(Genus)]
meta_top[, rank := .I]

if (nrow(meta_top) == 0) stop("No top genera selected.")

# Merge rank into forest table
forest_top <- merge(
  forest,
  meta_top,
  by = "Genus",
  all = FALSE
)

if (nrow(forest_top) == 0) stop("No matching genera in forest table.")

# Order by rank then cohort
setorder(forest_top, rank, Cohort)

# Remove rank column from export
forest_top[, rank := NULL]

dir.create("results/sex_paper/forest", showWarnings=FALSE, recursive=TRUE)

out_path <- "results/sex_paper/forest/forest_top10.csv"
fwrite(forest_top, out_path)

cat("Wrote:", out_path, "\nTop genera:\n")
print(meta_top$Genus)
