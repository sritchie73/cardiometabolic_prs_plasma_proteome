library(data.table)

args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]

bim <- fread(bim_file, header=FALSE)
bim[, varcount := seq_len(.N),by=.(V2)]
bim[, V2 := paste(V2, varcount, sep=":")] 
bim[, varcount := NULL]

fwrite(bim, file=bim_file, sep=" ",  quote=FALSE, col.names=FALSE)

