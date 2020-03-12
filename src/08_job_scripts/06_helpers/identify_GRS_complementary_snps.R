library(data.table)

args <- commandArgs(trailingOnly = TRUE)
complementary_snps_file <- args[1]
grs_file <- args[2]
out_file <- args[3]

complementary_snps <- fread(complementary_snps_file, header=FALSE, sep=" ") # otherwise detects ":" as the separator
setnames(complementary_snps, "rsid")

grs <- fread(grs_file)
cat(grs[complementary_snps, on = .(rsid), nomatch=0L, .N], "\n", sep="", file=out_file)
