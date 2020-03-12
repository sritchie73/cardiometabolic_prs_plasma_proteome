library(data.table)

args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]

# So we can handle multi-allelic variants, give them identifiers of the form
# rsID:A1:A2 where A1 and A2 are alphabetically sorted
bim <- fread(bim_file, header=FALSE)
setnames(bim, c("chr", "rsid", "cm", "pos", "A1", "A2"))
bim[, row := .I]
bim[rsid != ".", rsid := paste(rsid, paste(sort(c(A1, A2)), collapse=":"), sep=":"), by = row]
bim[, row := NULL]

# Create chr:pos:allele IDs for variants without rsIDs so plink doesn't complain
bim[rsid == ".", rsid := paste(chr, pos, A1, A2, sep=":")]

fwrite(bim, file=bim_file, sep="\t",  quote=FALSE, col.names=FALSE)
