library(data.table)

repo_root <- "C:/work/PD_microbiome_sex"
raw_dir   <- file.path(repo_root, "data/raw/malaysia_mothur/Malasia")
out_dir   <- file.path(repo_root, "data/processed/malaysia")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

clean_tax_token <- function(tok) {
  tok <- gsub('"', "", tok)
  tok <- gsub("\\(\\d+\\)", "", tok)
  trimws(tok)
}

extract_deepest_informative <- function(tax_string) {
  if (is.na(tax_string) || tax_string == "") return("Unclassified")

  parts <- unlist(strsplit(tax_string, ";", fixed = TRUE))
  parts <- gsub('"', "", parts)
  parts <- gsub("\\(\\d+\\)", "", parts)
  parts <- trimws(parts)
  parts <- parts[parts != ""]

  # walk from deepest to root
  for (tok in rev(parts)) {
    if (!grepl("unclassified", tok, ignore.case = TRUE) &&
        !grepl("^bacteria$", tok, ignore.case = TRUE)) {
      return(tok)
    }
  }

  return("Unclassified")
}

read_shared_mothur <- function(path) {
  x <- fread(path)
  if (!all(c("label","Group","numOtus") %in% names(x)))
    stop(paste("Unexpected .shared format in", path))
  otu_cols <- setdiff(names(x), c("label","Group","numOtus"))
  dt <- x[, c("Group", otu_cols), with = FALSE]
  setnames(dt, "Group", "SampleID")
  dt[, (otu_cols) := lapply(.SD, as.numeric), .SDcols = otu_cols]
  dt
}

# ---- READ COUNTS (PD vs Control) ----
pd_counts  <- read_shared_mothur(file.path(raw_dir, "FINAL_PD_TO.opti_mcc.shared"))
ctl_counts <- read_shared_mothur(file.path(raw_dir, "FINAL_CONTROL_TO.opti_mcc.shared"))
counts <- rbindlist(list(pd_counts, ctl_counts), use.names = TRUE, fill = TRUE)

otu_cols <- setdiff(names(counts), "SampleID")
counts[, (otu_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = otu_cols]

# PD label from membership
counts[, PD := fifelse(SampleID %in% pd_counts$SampleID, 1L, 0L)]

# Sex label from membership (male/female splits)
male_ids   <- read_shared_mothur(file.path(raw_dir, "FINAL_MALE_TO.opti_mcc.shared"))$SampleID
female_ids <- read_shared_mothur(file.path(raw_dir, "FINAL_FEMALE_TO.opti_mcc.shared"))$SampleID

counts[, Sex := fifelse(SampleID %in% male_ids, "Male",
                 fifelse(SampleID %in% female_ids, "Female", NA_character_))]

# Fail if missing Sex
if (any(is.na(counts$Sex))) {
  nmiss <- sum(is.na(counts$Sex))
  stop(paste("Malaysia: missing Sex for", nmiss, "samples. Check raw split files."))
}

# ---- TAXONOMY (use PD taxonomy file; should share OTU ids) ----
tax <- fread(file.path(raw_dir, "FINAL_PD_TO.opti_mcc.0.03.cons.taxonomy"))
tax[, Genus := vapply(Taxonomy, extract_deepest_informative, character(1))]

# ---- OTU -> GENUS collapse ----
genus_map <- tax[match(otu_cols, OTU), Genus]
genus_map[is.na(genus_map)] <- "Unclassified"

long <- melt(
  counts,
  id.vars = c("SampleID","PD","Sex"),
  measure.vars = otu_cols,
  variable.name = "OTU",
  value.name = "Count"
)

long[, Genus := genus_map[match(OTU, otu_cols)]]
genus_dt <- long[, .(Count = sum(Count, na.rm = TRUE)), by = .(SampleID, PD, Sex, Genus)]
wide <- dcast(genus_dt, SampleID + PD + Sex ~ Genus, value.var = "Count", fill = 0)

# Write outputs
fwrite(wide, file.path(out_dir, "malaysia_genus_counts_rebuilt.csv"))
fwrite(wide[, .(SampleID, PD, Sex)], file.path(out_dir, "malaysia_metadata_rebuilt.csv"))

cat("Malaysia rebuilt.\n")
cat("Samples:", nrow(wide), "\n")
cat("Genera:", ncol(wide) - 3, "\n")
print(wide[, .N, by=.(PD, Sex)][order(PD, Sex)])
