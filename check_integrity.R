library(data.table)

meta   <- fread("data/processed/ipd_sex/ipd_metadata.csv")
union  <- fread("data/processed/ipd_paper2/ipd_genus_counts_union_rebuilt.csv")
shared <- fread("data/processed/ipd_paper2/ipd_genus_counts_shared2_rebuilt.csv")

cat("Meta rows  :", nrow(meta), "\n")
cat("Union rows :", nrow(union), "\n")
cat("Shared rows:", nrow(shared), "\n")

stopifnot(
  all(meta$SampleID %in% union$SampleID),
  all(meta$SampleID %in% shared$SampleID),
  nrow(meta) == nrow(shared)
)

cat("Integrity OK\n")