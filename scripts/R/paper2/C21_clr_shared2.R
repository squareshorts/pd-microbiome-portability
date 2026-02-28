cat("C21 CLR STARTED\n")
cat("WD:", getwd(), "\n")

library(data.table)

# ---- Paths ----
in_path  <- "data/processed/ipd_paper2/ipd_genus_counts_shared2_rebuilt.csv"
out_path <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"

# ---- Load ----
dt <- fread(in_path)

cat("Input dim:", dim(dt), "\n")

# ---- Separate IDs ----
ids <- dt$SampleID
X   <- as.matrix(dt[, !"SampleID"])

mode(X) <- "numeric"

# ---- Sanity checks ----
if (any(is.na(X))) stop("NA values detected in counts matrix")

if (any(X < 0)) stop("Negative counts detected")

cat("Min count:", min(X), "\n")
cat("Max count:", max(X), "\n")

# ---- Add pseudocount ----
# Small constant to avoid log(0)
X <- X + 1

# ---- CLR transformation ----
# geometric mean per sample
gm <- exp(rowMeans(log(X)))

# divide by geometric mean
X_clr <- log(X / gm)

# ---- Reassemble ----
clr_dt <- data.table(SampleID = ids)
clr_dt <- cbind(clr_dt, as.data.table(X_clr))

# ---- Diagnostics ----
cat("CLR dim:", dim(clr_dt), "\n")

# Mean of each row should be ~0
row_means <- rowMeans(X_clr)
cat("Row mean range (should be near 0):\n")
print(range(row_means))

# ---- Save ----
fwrite(clr_dt, out_path)

cat("Saved to:", out_path, "\n")
cat("C21 CLR DONE\n")