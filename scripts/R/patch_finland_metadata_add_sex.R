library(data.table)

repo_root <- "C:/work/PD_microbiome_sex"
raw_dir   <- file.path(repo_root, "data/raw/finland_mothur/Finlandia")
proc_dir  <- file.path(repo_root, "data/processed/finland")

meta_path <- file.path(proc_dir, "finland_metadata.csv")
meta <- fread(meta_path)

# --- Identify male / female sample IDs from raw mothur splits ---
read_shared_ids <- function(path) {
  x <- fread(path)
  x$Group
}

male_ids <- read_shared_ids(file.path(raw_dir, "FINAL_FIN_male.opti_mcc.shared"))
female_ids <- read_shared_ids(file.path(raw_dir, "FINAL_FIN_FEMALE.opti_mcc.shared"))

meta[, Sex := fifelse(SampleID %in% male_ids, "Male",
               fifelse(SampleID %in% female_ids, "Female", NA_character_))]

# Sanity checks
cat("Total samples:", nrow(meta), "\n")
cat("Missing Sex:", sum(is.na(meta$Sex)), "\n")
print(meta[, .N, by = .(PD, Sex)][order(PD, Sex)])

# Fail if any NA
if (any(is.na(meta$Sex))) {
  stop("Some Finland samples could not be assigned Sex.")
}

# Write patched file
fwrite(meta, file.path(proc_dir, "finland_metadata_sex.csv"))
cat("Wrote finland_metadata_sex.csv\n")
