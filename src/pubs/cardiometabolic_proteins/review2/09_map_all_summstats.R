library(data.table)
library(openxlsx)
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

prs <- c("CAD_metaGRS", "Stroke_metaGRS", "Afib_2018", "CKD_2019", "T2D_2018")

# Load GWASs of interest
gwass <- sprintf("analyses/mendelian_randomisation/%s/gwas_summary_stats", prs)
gwass <- sapply(gwass, list.files, full.names=TRUE)
if(length(gwass) == 0) stop(sprintf("No gwas summary stats found in analyses/mendelian_randomisation/%s/gwas_summary_stats", prs))
gwass <- data.table(PRS=prs, gwas_dir=gwass, gwas_name=basename(gwass))

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

# Flag ambiguous SNPs for removal
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

# Remove aptamers with comparable binding to other proteins or isoforms
#
# N = 39 removed aptamers, 34 removed proteins
# N = 922 aptamers remaining, 856 proteins
# 
prot_info <- prot_info[!(
  tolower(Cross_Reactivity) %like% "comparable binding" | 
  tolower(Cross_Reactivity) %like% "similar affinity"
)]

# Get all pQTLs that can be mapped to GWAS summary stats
apt_ivs <- foreach(gwas_idx = gwass[,.I], .combine=rbind) %do% {
  cat(sprintf("Loading GWAS summary statistics for GWAS %s for PRS %s...\n", gwass[gwas_idx, gwas_name], gwass[gwas_idx, PRS]))
  # Load GWAS summary stats
  gwas_ss <- foreach(chr_idx = 1:22, .combine=rbind) %do% {
    chr_ss <- fread(sprintf("%s/chr%s.txt.gz", gwass[gwas_idx, gwas_dir], chr_idx))
    chr_ss <- chr_ss[pos %in% cis_pQTLs[chr == chr_idx, pos]]
    chr_ss
  }
  
  # Flag strand ambiguous alleles, drop if no MAF provided
  gwas_ss[, ambig := EA == flip_strand(OA)]
  if(!("MAF" %in% names(gwas_ss))) {
    gwas_ss[, MAF := NA] # placeholder value
  }
  gwas_ss <- gwas_ss[!(ambig) | (!is.na(MAF) & (MAF < ambig_max_maf | MAF > 1 - ambig_max_maf))]

	# Filter to biallelic SNPs if requested
	if (biallelic_only) {
		multi <- unique(gwas_ss[, .(chr, pos, EA, OA)])
		multi <- multi[,.N,by=.(chr, pos)][N > 1]
		gwas_ss <- gwas_ss[!multi, on = .(chr, pos)]
	}

	if (snps_only) {
		gwas_ss <- gwas_ss[nchar(EA) == 1 & nchar(OA) == 1]
	}

  cat("Matching GWAS to pQTL summary statistics...\n")
  
  # Map pQTLs and filter to independent instruments
  apt_ivs <- foreach(this_aptamer = prot_info$SOMAMER_ID, .combine=rbind) %dopar% {
    # Map cis-pQTLs to GWAS summary stats
    this_cis_pqtls <- cis_pQTLs[SOMAMER_ID == this_aptamer]
    cis_ivs <- this_cis_pqtls[gwas_ss, on = .(chr, pos), nomatch=0]
    setnames(cis_ivs, gsub("i\\.", "gwas.", names(cis_ivs)))
    if(nrow(cis_ivs) == 0) return(NULL)

    # Check for matching by alleles, flipping strand if necessary
    cis_ivs[, match := EA == gwas.EA | EA == gwas.OA]
    cis_ivs[!(match), c("gwas.EA", "gwas.OA") := .(flip_strand(gwas.EA), flip_strand(gwas.OA))]
    cis_ivs[!(match), match := EA == gwas.EA | EA == gwas.OA]
    cis_ivs <- cis_ivs[(match)]
    if(nrow(cis_ivs) == 0) return(NULL)

    # Orient effect alleles
    cis_ivs[EA != gwas.EA, c("gwas.EA", "gwas.OA", "MAF", "effect") := .(gwas.OA, gwas.EA, 1 - MAF, -effect)]
   
    # Make sure ambiguous alleles are correctly oriented, if MAF has been provided
		cis_ivs[((ambig) | (gwas.ambig)) & ((EAF < 0.5 & MAF > 0.5) | (EAF > 0.5 & MAF < 0.5)),
						c("gwas.EA", "gwas.OA", "MAF", "effect") := .(gwas.OA, gwas.EA, 1 - MAF, -effect)]
    if(nrow(cis_ivs) == 0) return(NULL)

    # Reorganise and rename columns
    cis_ivs[,.(PRS = gwass[gwas_idx, PRS], GWAS = gwass[gwas_idx, gwas_name], Aptamer = this_aptamer, 
							 pQTL.chr=chr, pQTL.pos=pos, pQTL.EA=EA, pQTL.OA=OA, pQTL.EAF=EAF,
							 pQTL.beta=effect, pQTL.SE=se, pQTL.P=P, gwas.beta=gwas.effect,
							 gwas.se, gwas.P, gwas.EAF=MAF)]
  }
}

# There are 916 / 922 aptamers that had cis-PQTls that could be mapped to GWAS summary stats
# covering 850 / 856 proteins

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
setkey(apt_ivs, pQTL.chr, pQTL.pos)

cat("Pruning cis-variants to those independent by LD...\n")

# Take the top variant (by effect and p-value) and add it to the list of independent IVs
# then remove from apt_ivs all remaining SNPs in LD, repeating until no SNPs are left
ind_ivs <- apt_ivs[0]
iter = 1
while(nrow(apt_ivs) > 0) {
  cat(sprintf("Iteration %s, IVs remaining: %s\n", iter, apt_ivs[,.N]))
  top_cis <- apt_ivs[order(-abs(pQTL.beta))][order(pQTL.P)][order(Aptamer)][order(GWAS)][order(PRS)][, .SD[1], by=.(PRS, GWAS, Aptamer)]
  ind_ivs <- rbind(ind_ivs, top_cis)

  top_cis <- top_cis[, .(PRS, GWAS, Aptamer, pQTL.chr, pQTL.pos)]
  setkey(top_cis, pQTL.chr, pQTL.pos)
  ld_with_top <- rbind(
    cis_ld[top_cis, on = .(chromosome=pQTL.chr, position1=pQTL.pos), nomatch=0][,.(PRS, GWAS, Aptamer, chr=chromosome, pos=position2)],
    cis_ld[top_cis, on = .(chromosome=pQTL.chr, position2=pQTL.pos), nomatch=0][,.(PRS, GWAS, Aptamer, chr=chromosome, pos=position1)]
  )
 
  apt_ivs <- apt_ivs[!ld_with_top, on = .(PRS, GWAS, Aptamer, pQTL.chr=chr, pQTL.pos=pos)]
  apt_ivs <- apt_ivs[!top_cis, on =.(PRS, GWAS, Aptamer, pQTL.chr, pQTL.pos)]
  iter <- iter + 1
}

# Filter to aptamers with at least 3 independent (by-LD) cis-pQTLs for MR:
# 
# 536 / 916 aptamers, covering 497 proteins
#
pass <- ind_ivs[,.N,by=.(PRS, GWAS, Aptamer)][N >= 3, .(PRS, GWAS, Aptamer)]
ind_ivs <- ind_ivs[pass, on = .(PRS, GWAS, Aptamer)]

# Orient effects to minor allele in INTERVAL
ind_ivs[pQTL.EAF > 0.5, c("pQTL.EA", "pQTL.OA", "pQTL.EAF", "pQTL.beta", "gwas.beta", "gwas.EAF") := 
                        .(pQTL.OA, pQTL.EA, 1 - pQTL.EAF, -pQTL.beta, -gwas.beta, 1-gwas.EAF)]

# Add protein information
setnames(ind_ivs, "Aptamer", "SOMAMER_ID")
ind_ivs <- prot_info[ind_ivs, on = .(SOMAMER_ID), nomatch=0]

# Select columns
ind_ivs <- ind_ivs[, .(PRS, GWAS, Target, UniProt, Gene, Gene.chr=chr, Gene.start=start, 
                       Aptamer, pQTL.chr, pQTL.pos, pQTL.EA, pQTL.OA, pQTL.EAF, 
                       pQTL.beta, pQTL.SE, pQTL.P,
                       gwas.EAF, gwas.beta, gwas.SE=gwas.se, gwas.P)]

# Order rows
ind_ivs <- ind_ivs[order(pQTL.P)][order(Aptamer)][order(Target)][order(GWAS)][order(PRS)]

# Write out
fwrite(ind_ivs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/cis_ivs.txt")


