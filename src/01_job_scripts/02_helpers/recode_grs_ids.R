library(data.table)

args <- commandArgs(trailingOnly = TRUE)
grs_file <- args[1]

# columns should be: rsid, chr, pos, effect_allele, other_allele, and weight
grs <- fread(grs_file)
expected <- c("rsid", "chr", "pos", "effect_allele", "other_allele", "weight")
if (length(union(names(grs), expected)) != length(expected)) {
  stop("Expected columns in GRS:\n", paste(paste("'", expected, "'"), collapse=", "),
       "\nFound:", paste(paste("'", names(grs), "'"), collapse=", "))
}

# Create chr:pos:allele ids. In this format alleles are sorted
# alphabetically. This is the most straighforward way to enable
# string matching of IDs in plink, while handling multi-allelic
# variants, without having to wrangle issues of which order the
# alleles appear in any given dataset.
grs[, row := .I]
grs[, sorted_alleles := paste(sort(c(effect_allele, other_allele)), collapse=":"), by=row]
grs[, new_id := paste(chr, pos, sorted_alleles, sep=":")]

grs <- grs[, .(rsid=new_id, chr, pos, effect_allele, other_allele, weight)]
out_fn <- gsub(".gz$", "", grs_file)
fwrite(grs, file=out_fn, sep="\t",  quote=FALSE)
system(sprintf("gzip -f %s", out_fn), wait=TRUE)
