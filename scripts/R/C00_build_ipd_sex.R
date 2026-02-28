library(data.table)

cohorts <- list(
  finland  = list(counts="data/processed/finland/finland_genus_counts.csv",
                  meta  ="data/processed/finland/finland_metadata_sex.csv"),
  malaysia = list(counts="data/processed/malaysia/malaysia_genus_counts.csv",
                  meta  ="data/processed/malaysia/malaysia_metadata.csv"),
  usa      = list(counts="data/processed/usa/usa_genus_counts.csv",
                  meta  ="data/processed/usa/usa_metadata.csv")
)

read_and_align <- function(cohort_name, paths) {

  X <- fread(paths$counts)
  M <- fread(paths$meta)

  # Standardize metadata names
  setnames(M, tolower(names(M)))

  if (!("sampleid" %in% names(M))) stop(paste(cohort_name, "metadata lacks SampleID"))
  if (!("pd" %in% names(M))) stop(paste(cohort_name, "metadata lacks PD"))
  if (!("sex" %in% names(M))) stop(paste(cohort_name, "metadata lacks Sex"))

  # Normalize sex
  M[, sex := fifelse(tolower(sex) %in% c("m","male"), "Male",
              fifelse(tolower(sex) %in% c("f","female"), "Female", sex))]

  # Normalize PD
  if (is.character(M$pd)) {
    M[, pd := fifelse(tolower(pd) %in% c("pd","parkinson","parkinsons","1","true"), 1L, 0L)]
  } else {
    M[, pd := as.integer(pd)]
  }

  M[, cohort := cohort_name]

  # Ensure counts table has SampleID
  if (!("SampleID" %in% names(X))) {
    nm <- names(X)
    idx <- which(tolower(nm)=="sampleid")
    if (length(idx)==1) setnames(X, nm[idx], "SampleID")
  }
  if (!("SampleID" %in% names(X))) stop(paste(cohort_name, "counts lacks SampleID"))

  # Join metadata onto counts
  setkey(X, SampleID)
  setkey(M, sampleid)

  X2 <- X[M, nomatch=0]

  # Remove join key duplication
  X2[, sampleid := NULL]

  # Build standardized outputs
  meta_out <- X2[, .(SampleID, PD=pd, Sex=sex, Cohort=cohort)]

  meta_cols <- c("SampleID","pd","sex","cohort","TYPE","type")
  genus_cols <- setdiff(names(X2), meta_cols)

  counts_out <- X2[, c("SampleID", genus_cols), with=FALSE]

  list(meta=meta_out, counts=counts_out)
}

parts <- lapply(names(cohorts), function(nm) read_and_align(nm, cohorts[[nm]]))

meta_all <- rbindlist(lapply(parts, `[[`, "meta"), use.names=TRUE, fill=TRUE)
counts_all <- rbindlist(lapply(parts, `[[`, "counts"), use.names=TRUE, fill=TRUE)

# Harmonize genus set using intersection
genus_sets <- lapply(parts, function(p) setdiff(names(p$counts), "SampleID"))
genus_common <- Reduce(intersect, genus_sets)

cat("Common genera:", length(genus_common), "\n")

counts_common <- counts_all[, c("SampleID", genus_common), with=FALSE]

dir.create("data/processed/ipd_sex", showWarnings = FALSE, recursive = TRUE)

fwrite(meta_all, "data/processed/ipd_sex/ipd_metadata.csv")
fwrite(counts_common, "data/processed/ipd_sex/ipd_genus_counts_common.csv")

cat("IPD samples:", nrow(meta_all), "\n")
print(meta_all[, .N, by=.(Cohort, PD, Sex)][order(Cohort, PD, Sex)])
