library(data.table)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
omic_dir <- args[1]

# combine per-chromosome dosages, also split out file to know
# which allele is the dosage allele later
files <- list.files(pattern="*.raw", path=omic_dir, full.names=TRUE)
geno_long <- rbindlist(lapply(files, function(x) { 
  dt <- fread(x)
  melt(dt, id.vars=names(dt)[1:6], variable.name="variant", value.name="dosage")
}))
geno_long <- geno_long[, .(IID, variant, dosage)]
dosage_alleles <- unique(geno_long[,.(variant)])
dosage_alleles[, allele := gsub(".*_", "", variant)]
dosage_alleles[, variant := gsub("_.*", "", variant)]
geno <- dcast(geno_long, IID ~ variant, value.var="dosage")
colnames(geno) <- gsub("_.*$", "", colnames(geno))
fwrite(geno, file=file.path(omic_dir, "qtl_geno_prob.tsv"), sep="\t", quote=FALSE)  
fwrite(dosage_alleles, file=file.path(omic_dir, "qtl_dosage_allele_prob.tsv"), sep="\t", quote=FALSE)

# Combine per-chromosome alleles
peds <- list.files(pattern="*\\.ped", path=omic_dir, full.names=TRUE)
maps <- list.files(pattern="*\\.map", path=omic_dir, full.names=TRUE)
merge_alleles <- function(x, y) { merge(x, y, by = "IID") }
alleles <- foreach(idx = seq_along(peds), .combine=merge_alleles) %do% {
  ped <- fread(peds[idx], header=FALSE)
  map <- fread(maps[idx], header=FALSE)
  setnames(map, c("chr", "rsid", "cm", "pos"))
  setnames(ped, c("FID", "IID", "Father", "Mother", "sex", "phenotype",
           paste(rep(map$rsid, each=2), rep(c("A1", "A2"), times=nrow(map)), sep="_")))
  ped <- melt(ped, id.vars=c("FID", "IID", "Father", "Mother", "sex", "phenotype"))
	ped[, rsid := gsub("_.*", "", variable)]
  ped[, allele := gsub(".*_", "", variable)]
  ped <- dcast(ped, IID + rsid ~ allele, value.var="value") 
  ped[, alleles := paste(A1, A2, sep="/")]
  # When heterozygous, make sure only one type (i.e. A/G and G/A become A/G).
  alleles <- unique(ped, by=c("rsid", "alleles"))
  hetero <- alleles[A1 != A2]
  hetero[, hetero_allele := paste(sort(c(A1,A2)), collapse="/"), by=1:nrow(hetero)]
  ped[hetero, on = .(rsid, alleles), alleles := hetero_allele]
  ped <- dcast(ped, IID ~ rsid, value.var="alleles")
  ped
} 
alleles <- alleles[geno[,.(IID)], on = .(IID)]
fwrite(alleles, file=file.path(omic_dir, "qtl_alleles.tsv"), sep="\t", quote=FALSE)

# also write out a list of variants that could not be matched / did not pass QC
varlist <- fread(sprintf("%s/qtl_list.txt", omic_dir))
geno_vars <- colnames(geno)[-1] #drop IID
missing_vars <- setdiff(varlist$variants, geno_vars)
if (length(missing_vars) > 0) {
  fwrite(data.table(rsID=missing_vars), file=sprintf("%s/missing_qtl_variants_prob.txt", omic_dir))
}
