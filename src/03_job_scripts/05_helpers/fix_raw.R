library(data.table)

args <- commandArgs(trailingOnly = TRUE)
raw_file <- args[1]
chr <- args[2]
pos <- args[3]

raw_dt <- fread(raw_file)
allele <- gsub("._", "", colnames(raw_dt)[7])
colnames(raw_dt)[7] <- paste0("chr", chr, ":", pos, "_", allele)

fwrite(raw_dt, sep="\t", quote=FALSE, file=raw_file)
