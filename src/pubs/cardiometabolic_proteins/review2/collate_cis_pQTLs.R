library(data.table)
library(foreach)
library(doMC)
source("src/utilities/flip_strand.R")

# Set up parallel environment
if (!exists("ncores")) {
  ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
  ncores <- as.integer(ncores)
  if(is.na(ncores)) ncores <- 1
}
registerDoMC(ncores)
setDTthreads(ncores)

# Ambiguous SNP MAF cutoff
ambig_max_maf <- 0.42 # set to 0 to drop ambiguous SNPs altogether

# Cis LD prune threshold
ld_thresh <- 0.1

# Other options
biallelic_only <- TRUE
snps_only <- TRUE

cat("Loading pQTL summary statistics...\n")

# Load cis-pQTLs
cis_pQTLs <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")

# Obtain their MAF
if (!("EAF" %in% names(cis_pQTLs))) {
  chr_stats <- foreach(chr_idx = 1:22, .combine=rbind) %do% {
		fread(sprintf("data/INTERVAL/reference_files/imputed_genotypes/impute_%s_interval.snpstats", chr_idx))
	}
  cis_pQTLs[chr_stats, on = .(chr=chromosome, pos=position, EA=minor_allele), EAF := MAF]
  cis_pQTLs[chr_stats, on = .(chr=chromosome, pos=position, EA=major_allele), EAF := 1-MAF]
  fwrite(cis_pQTLs, sep="\t", quote=FALSE, file="analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
}

# Remove ambiguous SNPs
cis_pQTLs[, ambig := EA == flip_strand(OA)]
cis_pQTLs <- cis_pQTLs[!(ambig) | EAF < ambig_max_maf | EAF > 1 - ambig_max_maf]

# Filter to biallelic SNPs if requested
if (biallelic_only) {
  multi <- unique(cis_pQTLs[, .(chr, pos, EA, OA)])
  multi <- multi[,.N,by=.(chr, pos)][N > 1]
  cis_pQTLs <- cis_pQTLs[!multi, on = .(chr, pos)]
}

if (snps_only) {
  cis_pQTLs <- cis_pQTLs[nchar(EA) == 1 & nchar(OA) == 1]
}

# Load protein information.
prot_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

# Filter to aptamers passing QC
prot_info <- prot_info[Type == "Protein"]

# Select columns
prot_info <- prot_info[, .(SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot,
                           Gene=Gene.Name, chr, start, Cross_Reactivity=Characterization.Info)]

# Fix bad entries (Aptamers for the same target with different/missing gene/uniprot information)
prot_info[Target == "14-3-3 protein family", UniProt := "P61981|Q04917"]
prot_info[Target == "Induced myeloid leukemia cell differentiation protein Mcl-1",
  c("UniProt", "Gene", "chr", "start") :=
  .("Q07820", "MCL1", "1", "150547027")]
prot_info[Target == "Protein delta homolog 1",
  c("UniProt", "Gene", "chr", "start") :=
  .("P80370", "DLK1", "14", "101193202")]
prot_info[Target == "Stromal cell-derived factor 1",
  c("UniProt", "Gene", "chr", "start") :=
  .("P48061", "CXCL12", "10", "44865601")]

# Filter to aptamers with cis-pQTLs:
prot_info <- prot_info[SOMAMER_ID %in% unique(cis_pQTLs$SOMAMER_ID)] # N = 961 aptamers, 890 proteins
cis_pQTLs <- cis_pQTLs[SOMAMER_ID %in% prot_info$SOMAMER_ID]

cat("Loading LD between cis-pQTLs for LD-pruning cis-IVs...\n")

# Load LD between cis-pQTLs so we can prune mapped cis-pQTLs for LD
if (!file.exists("analyses/mendelian_randomisation/pqtls/collated_cis_LD.txt.gz")) {
	cis_ld <- list.files("analyses/mendelian_randomisation/pqtls/SOMAMERs_cis_LD/", full.names=TRUE)
	cis_ld <- rbindlist(lapply(cis_ld, fread), fill=TRUE)
  cis_ld <- cis_ld[order(position2)][order(position1)][order(chromosome)]
	cis_ld <- unique(cis_ld)
  cis_ld[, r2 := correlation**2]
  fwrite(cis_ld, sep="\t", quote=FALSE, compress="gzip", file="analyses/mendelian_randomisation/pqtls/collated_cis_LD.txt.gz")
  fwrite(cis_ld[r2 >= 0.1], sep="\t", quote=FALSE, compress="gzip", file="analyses/mendelian_randomisation/pqtls/collated_cis_LD_r20.1.txt.gz")
} else {
  if (ld_thresh < 0.1) {
		cis_ld <- fread("analyses/mendelian_randomisation/pqtls/collated_cis_LD.txt.gz")
  } else {
    cis_ld <- fread("analyses/mendelian_randomisation/pqtls/collated_cis_LD_r20.1.txt.gz")
  }
}
cis_ld <- cis_ld[r2 > ld_thresh]
setkey(cis_ld, chromosome, position1, position2)

# Prepare apt_ivs for joining to cis_ld to filter cis_ivs by LD
setkey(cis_pQTLs, chr, pos)

cat("Pruning cis-variants to those independent by LD...\n")

# Take the top variant (by effect and p-value) and add it to the list of independent IVs
# then remove from apt_ivs all remaining SNPs in LD, repeating until no SNPs are left
ind_pQTLs <- cis_pQTLs[0]
iter = 1
while(nrow(cis_pQTLs) > 0) {
  cat(sprintf("Iteration %s, variants remaining: %s\n", iter, cis_pQTLs[,.N]))
  top_cis <- cis_pQTLs[order(-abs(effect))][order(P)][order(SOMAMER_ID)][, .SD[1], by=.(SOMAMER_ID)]
  ind_pQTLs <- rbind(ind_pQTLs, top_cis)

  top_cis <- top_cis[, .(SOMAMER_ID, chr, pos)]
  setkey(top_cis, chr, pos)
  ld_with_top <- rbind(
    cis_ld[top_cis, on = .(chromosome=chr, position1=pos), nomatch=0][,.(SOMAMER_ID, chr=chromosome, pos=position2)],
    cis_ld[top_cis, on = .(chromosome=chr, position2=pos), nomatch=0][,.(SOMAMER_ID, chr=chromosome, pos=position1)]
  )
 
  cis_pQTLs <- cis_pQTLs[!ld_with_top, on = .(SOMAMER_ID, chr, pos)]
  cis_pQTLs <- cis_pQTLs[!top_cis, on =.(SOMAMER_ID, chr, pos)]
  iter <- iter + 1
}

# Orient effects to minor allele in INTERVAL
ind_pQTLs[EAF > 0.5, c("EA", "OA", "EAF", "effect") := .(OA, EA, 1-EAF, -effect)]

# Add protein information
ind_pQTLs <- prot_info[ind_pQTLs, on = .(SOMAMER_ID), nomatch=0]

# Select columns
ind_pQTLs <- ind_pQTLs[, .(Target, UniProt, Gene, chr, start, Aptamer, pQTL.chr=i.chr, pQTL.pos=pos, EA,
                           OA, EAF, beta=effect, se, P, local_bonf, global_fdr)]

# Order rows
ind_pQTLs <- ind_pQTLs[order(P)][order(pQTL.chr)][order(Aptamer)][order(Target)]

# Write out
fwrite(ind_pQTLs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/cis_pQTLs.txt")


