cat("C20b CLEAN STARTED\n")
cat("SCRIPT WD:", getwd(), "\n")
library(data.table)

meta_path <- "data/processed/ipd_sex/ipd_metadata.csv"
uni_path  <- "data/processed/ipd_paper2/ipd_genus_counts_union_rebuilt.csv"

meta <- fread(meta_path)
Xuni <- fread(uni_path)

cat("META PATH:", normalizePath(meta_path), "\n")
cat("First 5 META IDs:\n")
print(head(meta$SampleID, 5))

cat("UNION PATH:", normalizePath(uni_path), "\n")
cat("First 5 UNION IDs:\n")
print(head(Xuni$SampleID, 5))

# Explicit value-based merge (no keys)
dt <- merge(meta[, .(SampleID)], Xuni, 
            by = "SampleID", 
            all.x = TRUE, 
            sort = FALSE)

cat("First 5 DT SampleID values:\n")
print(head(dt$SampleID, 5))

stopifnot(nrow(dt) == nrow(meta))

genus_cols <- setdiff(names(dt), "SampleID")

cohorts <- sort(unique(meta$Cohort))

present_in_cohort <- function(coh) {
  ids <- meta[Cohort == coh, SampleID]
  d <- dt[dt$SampleID %in% ids, ]
  colSums(as.matrix(d[, genus_cols, with = FALSE]) > 0) > 0
}

pres_mat <- sapply(cohorts, present_in_cohort)
n_coh_present <- rowSums(pres_mat)

keep <- names(n_coh_present)[n_coh_present >= 2]

cat("Union genera:", length(genus_cols), "\n")
cat("Keep (>=2 cohorts):", length(keep), "\n")

Xshared2 <- dt[, c("SampleID", keep), with = FALSE]

cat("First 5 IDs in Xshared2:\n")
print(head(Xshared2$SampleID, 5))

out_path <- "data/processed/ipd_paper2/ipd_genus_counts_shared2_rebuilt.csv"
fwrite(Xshared2, out_path)

cat("Dim:", dim(Xshared2), "\n")
cat("C20b CLEAN DONE\n")