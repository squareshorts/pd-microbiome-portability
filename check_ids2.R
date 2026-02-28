library(data.table)

cat("Working dir:", getwd(), "\n")

cat("Counts exists:", file.exists("data/processed/ipd_paper2/ipd_genus_counts_shared2_rebuilt.csv"), "\n")
cat("Meta exists  :", file.exists("data/processed/ipd_sex/ipd_metadata.csv"), "\n")

counts <- fread("data/processed/ipd_paper2/ipd_genus_counts_shared2_rebuilt.csv")
meta   <- fread("data/processed/ipd_sex/ipd_metadata.csv")

cat("Counts rows =", nrow(counts), "\n")
cat("Meta rows   =", nrow(meta), "\n")
cat("Intersection=", length(intersect(counts, meta)), "\n")
