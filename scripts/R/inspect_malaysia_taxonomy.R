library(data.table)

tax <- fread("data/raw/malaysia_mothur/Malasia/FINAL_PD_TO.opti_mcc.0.03.cons.taxonomy")

cat("Columns:\n")
print(names(tax))

cat("\nExample taxonomy strings:\n")
print(tax[1:20, .(OTU, Taxonomy)])

cat("\nLast-token frequency (raw):\n")
last_tok <- vapply(strsplit(tax$Taxonomy, ";", fixed=TRUE),
                   function(x) tail(x, 1), character(1))
last_tok <- gsub("\\(\\d+\\)", "", last_tok)
last_tok <- trimws(last_tok)
print(head(sort(table(last_tok), decreasing=TRUE), 30))
