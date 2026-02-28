library(data.table)

f <- fread("data/processed/finland/finland_genus_counts_rebuilt.csv")
m <- fread("data/processed/malaysia/malaysia_genus_counts_rebuilt.csv")
u <- fread("data/processed/usa/usa_genus_counts_rebuilt.csv")

gf <- setdiff(names(f), "SampleID")
gm <- setdiff(names(m), "SampleID")
gu <- setdiff(names(u), "SampleID")

cat("Finland genera:", length(gf), "\n")
cat("Malaysia genera:", length(gm), "\n")
cat("USA genera:", length(gu), "\n")

cat("Finland ∩ Malaysia:", length(intersect(gf, gm)), "\n")
cat("Finland ∩ USA:", length(intersect(gf, gu)), "\n")
cat("Malaysia ∩ USA:", length(intersect(gm, gu)), "\n")
cat("All three:", length(Reduce(intersect, list(gf, gm, gu))), "\n")
