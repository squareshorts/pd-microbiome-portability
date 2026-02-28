library(data.table)

cohorts <- list(
  finland  = "data/processed/finland/finland_genus_counts_rebuilt.csv",
  malaysia = "data/processed/malaysia/malaysia_genus_counts_rebuilt.csv",
  usa      = "data/processed/usa/usa_genus_counts_rebuilt.csv"
)

parts <- lapply(names(cohorts), function(nm) {
  dt <- fread(cohorts[[nm]])
  dt[, Cohort := nm]
  dt
})

meta_all <- rbindlist(lapply(parts, function(dt)
  dt[, .(SampleID, PD, Sex, Cohort)]))

counts_all <- rbindlist(lapply(parts, function(dt)
  dt[, !c("PD","Sex","Cohort"), with=FALSE]), fill=TRUE)

# union of genera
genus_cols <- setdiff(names(counts_all), "SampleID")
counts_all[is.na(counts_all)] <- 0

dir.create("data/processed/ipd_sex", showWarnings = FALSE, recursive = TRUE)

fwrite(meta_all, "data/processed/ipd_sex/ipd_metadata.csv")
fwrite(counts_all, "data/processed/ipd_sex/ipd_genus_counts_union.csv")

cat("IPD samples:", nrow(meta_all), "\n")
cat("Total genera (union):", length(genus_cols), "\n")

print(meta_all[, .N, by=.(Cohort, PD, Sex)][order(Cohort, PD, Sex)])
