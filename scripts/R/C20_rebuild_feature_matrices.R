cat("C20 FINAL CLEAN STARTED\n")
library(data.table)

meta_path <- "data/processed/ipd_sex/ipd_metadata.csv"
fin_path  <- "data/processed/finland/finland_genus_counts_rebuilt.csv"
mal_path  <- "data/processed/malaysia/malaysia_genus_counts_rebuilt.csv"
usa_path  <- "data/processed/usa/usa_genus_counts_rebuilt.csv"

meta <- fread(meta_path)
fin  <- fread(fin_path)
mal  <- fread(mal_path)
usa  <- fread(usa_path)

clean_counts <- function(dt) {
  drop_exact <- c("PD","Sex","Cohort")
  drop_prefix <- "^i\\."
  
  keep <- names(dt)[
    !(names(dt) %in% drop_exact) &
    !grepl(drop_prefix, names(dt))
  ]
  
  dt[, ..keep]
}

fin <- clean_counts(fin)
mal <- clean_counts(mal)
usa <- clean_counts(usa)

stopifnot("SampleID" %in% names(fin))
stopifnot("SampleID" %in% names(mal))
stopifnot("SampleID" %in% names(usa))

fin_g <- setdiff(names(fin), "SampleID")
mal_g <- setdiff(names(mal), "SampleID")
usa_g <- setdiff(names(usa), "SampleID")

union_g  <- Reduce(union, list(fin_g, mal_g, usa_g))
common_g <- Reduce(intersect, list(fin_g, mal_g, usa_g))

cat("Fin genera:", length(fin_g), "\n")
cat("Mal genera:", length(mal_g), "\n")
cat("USA genera:", length(usa_g), "\n")
cat("Union genera:", length(union_g), "\n")
cat("Common genera:", length(common_g), "\n")

add_missing <- function(dt, genes) {
  miss <- setdiff(genes, names(dt))
  if (length(miss)) dt[, (miss) := 0]
  cols <- c("SampleID", genes)
  dt[, ..cols]
}

fin_u <- add_missing(fin, union_g)
mal_u <- add_missing(mal, union_g)
usa_u <- add_missing(usa, union_g)

ipd_union <- rbindlist(list(fin_u, mal_u, usa_u), use.names=TRUE, fill=TRUE)

ipd_union <- ipd_union[SampleID %in% meta$SampleID]

dir.create("data/processed/ipd_paper2", showWarnings=FALSE, recursive=TRUE)

cat("WRITING UNION NOW...\n")
fwrite(ipd_union,
       "data/processed/ipd_paper2/ipd_genus_counts_union_rebuilt.csv")
cat("WRITE COMPLETE\n")

cat("Union dim:", dim(ipd_union), "\n")
cat("C20 FINAL CLEAN DONE\n")
