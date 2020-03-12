library(data.table)

args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]
out_file <- args[2]

bim <- fread(bim_file, header=FALSE)
setnames(bim, c("chr", "id", "nosex", "pos", "A1", "A2"))

# Flip allele strands
flip_strand <- function(alleles) {
  alleles <- gsub("A", "V", alleles)
  alleles <- gsub("T", "X", alleles)
  alleles <- gsub("C", "Y", alleles)
  alleles <- gsub("G", "Z", alleles)
  alleles <- gsub("Z", "C", alleles)
  alleles <- gsub("Y", "G", alleles)
  alleles <- gsub("X", "A", alleles)
  alleles <- gsub("V", "T", alleles)
  return(alleles)
}

complementary <- bim[A1 == flip_strand(A2)]

fwrite(complementary[,.(id)], file=out_file, quote=FALSE, col.names=FALSE)

