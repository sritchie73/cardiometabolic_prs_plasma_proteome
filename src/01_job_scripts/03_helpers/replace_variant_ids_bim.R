library(data.table)

args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]

bim <- fread(bim_file, header=FALSE)
setnames(bim, c("chr", "rsID", "cm", "pos", "A1", "A2"))

# Create chr:pos:allele ids. In this format alleles are sorted
# alphabetically. This is the most straighforward way to enable
# string matching of IDs in plink, while handling multi-allelic 
# variants, without having to wrangle issues of which order the 
# alleles appear in any given dataset.
bim[, row := .I]
bim[, sorted_alleles := paste(sort(c(A1,A2)), collapse=":"), by=row]
bim[, new_id := paste(chr, pos, sorted_alleles, sep=":")]

bim <- bim[, .(chr, new_id, cm, pos, A1, A2)]
fwrite(bim, file=bim_file, sep="\t",  quote=FALSE, col.names=FALSE)

