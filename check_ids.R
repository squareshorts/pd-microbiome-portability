library(data.table)

counts <- fread("data/processed/ipd_paper2/ipd_genus_counts_shared2_rebuilt.csv")
meta   <- fread("data/processed/ipd_sex/ipd_metadata.csv")

counts <- as.character(counts)
meta   <- as.character(meta)

cat("Counts rows =", nrow(counts), "\n")
cat("Meta rows   =", nrow(meta), "\n")
cat("Intersection=", length(intersect(counts, meta)), "\n")
