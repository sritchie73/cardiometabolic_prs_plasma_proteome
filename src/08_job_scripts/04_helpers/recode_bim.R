library(data.table)

args <- commandArgs(trailingOnly=TRUE)
bim_file <- args[1]

bim <- fread(bim_file)
bim[, V2 := gsub(":[A-Z]*:[A-Z]*$", "", V2)]
fwrite(bim, sep=" ", quote=FALSE, col.names=FALSE, file=bim_file)

