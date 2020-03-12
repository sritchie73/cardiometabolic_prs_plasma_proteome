library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

bim <- fread(input_file, header=FALSE)
setnames(bim, c("chr", "varid", "sex", "pos", "a1", "a2"))
bim[varid == ".", varid := paste0("chr", chr, ":", pos)]
fwrite(bim, file=output_file, sep=" ", col.names=FALSE, quote=FALSE)
