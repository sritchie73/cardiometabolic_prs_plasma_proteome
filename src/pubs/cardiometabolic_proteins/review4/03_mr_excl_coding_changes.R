library(data.table)
library(httr)
library(jsonlite)
library(xml2)
library(foreach)
library(doMC)
library(MendelianRandomization)
library(coloc)
library(ggplot2)
source("src/utilities/prot_pval.R")
source("src/utilities/flip_strand.R")
source("src/07_job_scripts/07_helpers/mr_functions.R")

# Set up parallel environment
if (!exists("ncores")) {
  ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
  ncores <- as.integer(ncores)
  if(is.na(ncores)) ncores <- 1
}
registerDoMC(ncores)
setDTthreads(ncores)

# Redo the MR excluding SNPs in LD with protein coding changes which may indicate aptamer
# binding effects

# Load data
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
mr_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/collated_mr.txt")

# Extract protein level associations
prs_assocs <- unique(prs_assocs[,.(PRS, Gene, PRS.FDR=FDR)])
mr_assocs <- unique(mr_assocs[,.(GWAS, Gene, MR=exp(mr_estimate), MR.L95=exp(mr_L95), MR.U95=exp(mr_U95), MR.P=mr_pval,
                                 MR.pleiotropy=pleiotropy_pval, colocalization)])

# Map between outcomes
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  GWAS = c("Afib", "CAD", "CKD", "StrokeIS", "T2DadjBMI")
)

prs_assocs[map, on = .(PRS), GWAS := i.GWAS]
mr_assocs[map, on = .(GWAS), PRS := i.PRS]

# Filter to FDR < 0.05 PRS to protein associations
prs_assocs <- prs_assocs[PRS.FDR < 0.05, .(PRS, GWAS, Gene)]
mr_assocs <- mr_assocs[prs_assocs, on = .(PRS, GWAS, Gene), nomatch=0]

# Load cis-pQTLs for these proteins so we can check for aptamer binding effects using VEP
cis <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
prot_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
cis[prot_info, on = .(SOMAMER_ID), Gene := Gene.Name]
cis <- cis[Gene %in% prs_assocs$Gene]
cis_info <- prot_info[Gene.Name %in% mr_assocs$Gene]

# Get rsIDs
pvar <- foreach(this_chr = unique(cis$chr), .combine=rbind) %do% {
  this_pvar <- fread(sprintf("data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/plink_format/pgen/impute_dedup_%s_interval.pvar", this_chr))
}
cis[pvar, on = .(chr=`#CHROM`, pos=POS), rsID := ID]
cis <- cis[order(pos)][order(chr)]

# Look up cis-pQTLs in ENSEMBL variant effect predictor (VEP) on GRCh37 to identify
# aptamers whose cis-pQTLs suggest potential aptamer binding effects. Maximum size of 
# POST request is 300 variants

unlist_json <- function(l) {
  foreach(li = seq_along(l), .combine=c) %do% {
    if (is.null(l[[li]])) {
      return(NA)
    } else if (!is.null(dim(l[[li]]))) {
      return(paste(as.vector(l[[li]]), collapse=";"))
    } else {
      return(l[[li]])
    }
  }
}

rbindf <- function(...) { rbind(..., fill=TRUE) }

cis[, BATCH := floor(.I/300)+1]
vep <- foreach(post_batch = unique(cis$BATCH), .combine=rbind) %dopar% {
	r <- POST("https://grch37.rest.ensembl.org/vep/human/id", content_type("application/json"), accept("application/json"),
						body = sprintf('{ "ids" : [ %s ] }', paste(cis[BATCH == post_batch & rsID %like% "^rs", sprintf('"%s"', unique(rsID))], collapse=", ")))

  stop_for_status(r)

  # Extract JSON to data.table
  vep <- foreach(ii = seq_along(content(r)), .combine=rbindf) %do% {
		rsID <- content(r)[[ii]][["id"]]
		most_severe <- content(r)[[ii]][["most_severe_consequence"]]
		consequences <- content(r)[[ii]][["transcript_consequences"]]
		consequences <- fromJSON(toJSON(consequences))
		consequences <- as.data.table(lapply(consequences, unlist_json))
		cbind(ID = rsID, most_severe = most_severe, consequences)
  }

  # Filter to most severe consequences per ID
  vep <- vep[, .SD[consequence_terms %like% unique(most_severe)], by=ID]

  # Filter columns
  vep[, .(rsID=ID, Gene=gene_symbol, impact, most_severe)]
}
cis[vep, on = .(rsID, Gene), c("VEP_most_severe_consequence", "impact") := .(most_severe, impact)]

# Manually use https://grch37.ensembl.org/Homo_sapiens/Tools/VEP for handful of
# variants without rsIDs
cis[rsID == "1:196690953:C:CTC" & Gene == "CFH", VEP_most_severe_consequence := "intron_variant"]
cis[rsID == "15:80440838:C:T" & Gene == "FAH", VEP_most_severe_consequence := "upstream_gene_variant"] 
cis[rsID == "5:42697657:C:CA" & Gene == "GHR", VEP_most_severe_consequence := "intron_variant"]
cis[rsID == "5:42720068:T:TA" & Gene == "GHR", VEP_most_severe_consequence := "3_prime_UTR_variant"]

# Identify coding changes that may influence aptamer binding affinity:
# Those listed are the ones that occur in the set of Genes examined, for
# full list see:
# https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
coding <- cis[VEP_most_severe_consequence %in% c("missense_variant", "splice_acceptor_variant", "splice_region_variant")]

# Write out:
fwrite(coding, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/coding_mutations.txt")

# Load cis variant LD
cis_ld <- fread("analyses/mendelian_randomisation/pqtls/collated_cis_LD_r20.1.txt.gz")

# Get all variants with r2 > 0.1 for each coding variant
coding_ld <- unique(coding[,.(Gene, chr, pos, type=VEP_most_severe_consequence)])
coding_ld <- rbind(
  coding_ld[cis_ld, on = .(chr=chromosome, pos=position1), nomatch=0, .(Gene, chr, pos, type, pos_in_ld=position2, r2)],
  coding_ld[cis_ld, on = .(chr=chromosome, pos=position2), nomatch=0, .(Gene, chr, pos, type, pos_in_ld=position1, r2)]
)
coding_ld <- unique(coding_ld)

# Drop from candidate cis-pQTL table
cis <- cis[!coding, on = .(Gene, chr, pos)]
cis <- cis[!coding_ld, on = .(Gene, chr, pos=pos_in_ld)]

# Load GWASs of interest
prs <- c("CAD_metaGRS", "Stroke_metaGRS", "CKD_2019", "T2D_2018")
gwass <- sprintf("analyses/mendelian_randomisation/%s/gwas_summary_stats", prs)
gwass <- sapply(gwass, list.files, full.names=TRUE)
if(length(gwass) == 0) stop(sprintf("No gwas summary stats found in analyses/mendelian_randomisation/%s/gwas_summary_stats", prs))
gwass <- data.table(PRS=prs, gwas_dir=gwass, gwas_name=basename(gwass))

# Ambiguous SNP MAF cutoff
ambig_max_maf <- 0.42 # set to 0 to drop ambiguous SNPs altogether

# Other options
biallelic_only <- TRUE
snps_only <- TRUE

# Flag ambiguous SNPs for removal
cis[, ambig := EA == flip_strand(OA)]
cis <- cis[!(ambig) | EAF < ambig_max_maf | EAF > 1 - ambig_max_maf]

# Filter to biallelic SNPs 
if (biallelic_only) {
  multi <- unique(cis[, .(chr, pos, EA, OA)])
  multi <- multi[,.N,by=.(chr, pos)][N > 1]
  cis <- cis[!multi, on = .(chr, pos)]
}

if (snps_only) {
  cis <- cis[nchar(EA) == 1 & nchar(OA) == 1]
}

# Get MR tests for mapping
ivs <- fread("analyses/pub/cardiometabolic_proteins/review3/collated_ivs.txt")
tests <- unique(ivs[,.(PRS, GWAS, Gene)])
gwass[, PRS := fcase(PRS == "CAD_metaGRS", "CAD_PRS", PRS == "Stroke_metaGRS", "IS_PRS", 
                     PRS == "CKD_2019", "CKD_PRS", PRS == "T2D_2018", "T2D_PRS")]

# Get all pQTLs that can be mapped to GWAS summary stats
apt_ivs <- foreach(gwas_idx = gwass[,.I], .combine=rbind) %do% {
  # get all genes we're testing in MR for this GWAS
  this_genes <- tests[GWAS == gwass[gwas_idx, gwas_name]]

  # Get all cis pQTLs for the gene(s)
  this_cis <- cis[Gene %in% this_genes$Gene]

  # Load GWAS summary stats that map to cis-pQTL positions
  gwas_ss <- foreach(chr_idx = 1:22, .combine=rbind) %do% {
    chr_ss <- fread(sprintf("%s/chr%s.txt.gz", gwass[gwas_idx, gwas_dir], chr_idx))
    chr_ss <- chr_ss[pos %in% this_cis[chr == chr_idx, pos]]
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

  # Map pQTLs and filter to candidate cis-pQTLs instruments
  foreach(this_aptamer = cis_info[Gene.Name %in% this_genes$Gene, SOMAMER_ID], .combine=rbind) %dopar% {
    # Map cis-pQTLs to GWAS summary stats
    this_cis_pqtls <- cis[SOMAMER_ID == this_aptamer]
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
    cis_ivs[,.(PRS = gwass[gwas_idx, PRS], GWAS = gwass[gwas_idx, gwas_name], 
               Gene = cis_info[SOMAMER_ID == this_aptamer, Gene.Name], Aptamer = this_aptamer,
               pQTL.chr=chr, pQTL.pos=pos, pQTL.EA=EA, pQTL.OA=OA, pQTL.EAF=EAF,
               pQTL.beta=effect, pQTL.SE=se, pQTL.P=P, gwas.beta=gwas.effect,
               gwas.se, gwas.P, gwas.EAF=MAF)]
  }
}

# Prepare apt_ivs for joining to cis_ld to filter cis_ivs by LD
setkey(apt_ivs, pQTL.chr, pQTL.pos)

# Take the top variant (by effect and p-value) and add it to the list of independent IVs
# then remove from apt_ivs all remaining SNPs in LD, repeating until no SNPs are left
ind_ivs <- apt_ivs[0]
while(nrow(apt_ivs) > 0) {
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
}

# Filter to aptamers with at least 3 independent (by-LD) cis-pQTLs for MR:
pass <- ind_ivs[,.N,by=.(PRS, GWAS, Aptamer)][N >= 3, .(PRS, GWAS, Aptamer)]
ind_ivs <- ind_ivs[pass, on = .(PRS, GWAS, Aptamer)]

# Orient effects to minor allele in INTERVAL
ind_ivs[pQTL.EAF > 0.5, c("pQTL.EA", "pQTL.OA", "pQTL.EAF", "pQTL.beta", "gwas.beta", "gwas.EAF") :=
                        .(pQTL.OA, pQTL.EA, 1 - pQTL.EAF, -pQTL.beta, -gwas.beta, 1-gwas.EAF)]

# Add cis_info
ind_ivs <- ind_ivs[cis_info, on = .(Aptamer=SOMAMER_ID), nomatch=0]


# Select columns
ind_ivs <- ind_ivs[, .(PRS, GWAS, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, 
                       Gene, Gene.chr=chr, Gene.start=start,
                       Aptamer, pQTL.chr, pQTL.pos, pQTL.EA, pQTL.OA, pQTL.EAF,
                       pQTL.beta, pQTL.SE, pQTL.P,
                       gwas.EAF, gwas.beta, gwas.SE=gwas.se, gwas.P)]

# Order rows
ind_ivs <- ind_ivs[order(pQTL.P)][order(Aptamer)][order(Target)][order(GWAS)][order(PRS)]

# Write out
fwrite(ind_ivs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/cis_ivs.txt")

# Function to extract relevant coefficients from MR output
extract_coef <- function(mro, egger=FALSE) {
  if (!egger) {
    data.table(mr_estimate = mro@Estimate,
               mr_se = mro@StdError,
               mr_L95 = mro@CILower,
               mr_U95 = mro@CIUpper,
               mr_pval = mro@Pvalue)
  } else {
    data.table(mr_estimate = c(mro@Estimate, mro@Intercept),
               mr_se = c(mro@StdError.Int),
               mr_L95 = c(mro@CILower.Est, mro@CILower.Int),
               mr_U95 = c(mro@CIUpper.Est, mro@CIUpper.Int),
               mr_pval = c(mro@Causal.pval, mro@Pleio.pval))
  }
}

# Run MR
tests <- unique(ind_ivs[,.(PRS, GWAS, Target, UniProt, Gene, Aptamer)])
apt_mr <- foreach(test_idx = tests[,.I], .combine=rbind) %dopar% {
  # Get instruments for this aptamer and GWAS
  this_test <- tests[test_idx]
  this_ivs <- ind_ivs[this_test, on = .(PRS, GWAS, Target, UniProt, Gene, Aptamer)]

  # Construct MR Input object
  mri <- this_ivs[, mr_input(bx=pQTL.beta, bxse=pQTL.SE, by=gwas.beta, byse=gwas.SE)]

	# Run 5 MR methods
	ivw <- mr_ivw(mri)
	simple_median <- mr_median(mri, weighting="simple")
	weighted_median <- mr_median(mri, weighting="weighted")
	weighted_mode <- mr_mbe(mri, weighting="weighted", stderror="simple")
	egger <- mr_egger(mri)

	mro <- rbind(idcol="Method",
		"IVW"=extract_coef(ivw),
		"Median"=extract_coef(simple_median),
		"Weighted Median"=extract_coef(weighted_median),
		"Weighted Mode"=extract_coef(weighted_mode),
		"MR Egger"=extract_coef(egger, TRUE)
	)
	mro[6, Method := "(Intercept)"]

  mro <- cbind(this_test, mro)
  mro
}

# Average across aptamers
prot_mr <- apt_mr[, .(mr_estimate = mean(mr_estimate), mr_se = mean(mr_se),
                      mr_L95 = mean(mr_L95), mr_U95 = mean(mr_U95),
                      mr_pval = prot_pvalue(mr_pval, mr_estimate)),
                  by = .(PRS, GWAS, Target, UniProt, Gene, Method)]

# Take median of causal estimates
mr_summary <- prot_mr[Method != "(Intercept)",
    .(mr_estimate = median(mr_estimate), mr_se = median(mr_se),
      mr_L95 = median(mr_L95), mr_U95 = median(mr_U95),
      mr_pval = prot_pvalue(mr_pval, mr_estimate, median)),
  by = .(PRS, GWAS, Target, UniProt, Gene)]

# Row-order
mr_summary <- mr_summary[order(mr_pval)][order(GWAS)][order(PRS)]

# Add pleiotropy p-value
mr_summary[prot_mr[Method == "(Intercept)"], on = .(PRS, GWAS, Target, UniProt, Gene), pleiotropy_pval := i.mr_pval]

# Write out
fwrite(apt_mr, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/mr_all_aptamers.txt")
fwrite(prot_mr, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/mr_all_proteins.txt")
fwrite(mr_summary, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/mr_summary_all_proteins.txt")

# Colocalization analysis
coloc <- foreach(gwas_idx = gwass[,.I], .combine=rbind) %do% {
  this_prs <- gwass[gwas_idx, PRS]
  this_gwas <- gwass[gwas_idx, gwas_name]

  # Determine windows we need to load for this GWAS (200Kb around each pQTL)
  windows <- ind_ivs[PRS == this_prs & GWAS == this_gwas, .(
    chr = pQTL.chr, window_start = pQTL.pos - 200000, window_end = pQTL.pos + 200000
  )]

  # Load GWAS summary stats and filter to variants in the 200Kb windows
  gwas_ss <- foreach(chr_idx = unique(windows$chr), .combine=rbind) %dopar% {
    chr_ss <- fread(sprintf("%s/chr%s.txt.gz", gwass[gwas_idx, gwas_dir], chr_idx))
    chr_ss <- chr_ss[windows, on = .(chr, pos >= window_start, pos <= window_end), nomatch=0,
                     .(chr, pos=x.pos, EA, OA, MAF, effect, se, P)]
    unique(chr_ss) # windows may overlap
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

  # Flag ambiguous alleles
  gwas_ss[, ambig := EA == flip_strand(OA)]

  # Load INTERVAL variant statistics so we can get MAF
  int_maf <- foreach(chr_idx = unique(windows$chr), .combine=rbind) %dopar% {
    chr_ss <- fread(sprintf("data/INTERVAL/reference_files/imputed_genotypes/impute_%s_interval.snpstats", chr_idx))
    chr_ss <- chr_ss[,.(chr=chromosome, pos=position, minor_allele, MAF)]
    chr_ss <- chr_ss[windows, on = .(chr, pos >= window_start, pos <= window_end), nomatch=0,
                     .(chr, pos=x.pos, minor_allele, MAF)]
    unique(chr_ss) # windows may overlap
  }

  # Load information about GWAS
  gwas_info <- fread(sprintf("%s/info.txt", gwass[gwas_idx, gwas_dir]))

  # Run coloc for each pQTL
  this_cis_ivs <- ind_ivs[PRS == this_prs & GWAS == this_gwas]
  foreach(test_idx = this_cis_ivs[,.I], .combine=rbind) %dopar% {
    this_chr <- this_cis_ivs[test_idx, pQTL.chr]
    this_pos <- this_cis_ivs[test_idx, pQTL.pos]
    window_start <- this_pos - 200000
    window_end <- this_pos + 200000
    this_somaid <- this_cis_ivs[test_idx, Aptamer]

    # Load pQTL summary statistics
    pqtl_ss <- fread(sprintf("data/full_pQTL_summary_stats/%s/%s_chrom_%s_meta_1.tbl.gz", this_somaid, this_somaid, this_chr))
    pqtl_ss <- pqtl_ss[, .(chr=chromosome, pos=position, EA=toupper(Allele1), OA=toupper(Allele2), effect=Effect, se=StdErr, P=10^(`log(P)`))]
    pqtl_ss <- pqtl_ss[pos >= window_start & pos <= window_end]

    # Add MAF
    pqtl_ss[int_maf, on = .(chr, pos), EAF := ifelse(EA == minor_allele, MAF, 1 - MAF)]

    # Filter to biallelic SNPs if requested
    if (biallelic_only) {
      multi <- unique(pqtl_ss[, .(chr, pos, EA, OA)])
      multi <- multi[,.N,by=.(chr, pos)][N > 1]
      pqtl_ss <- pqtl_ss[!multi, on = .(chr, pos)]
    }

    if (snps_only) {
      pqtl_ss <- pqtl_ss[nchar(EA) == 1 & nchar(OA) == 1]
    }

    # Flag ambiguous alleles
    pqtl_ss[, ambig := EA == flip_strand(OA)]

    # Map to GWAS summary stats
    map_ss <- pqtl_ss[gwas_ss, on = .(chr, pos), nomatch=0]
    setnames(map_ss, gsub("i\\.", "gwas.", names(map_ss)))
    if(nrow(map_ss) == 0) return(NULL)

    # Check for matching by alleles, flipping strand if necessary
    map_ss[, match := EA == gwas.EA | EA == gwas.OA]
    map_ss[!(match), c("gwas.EA", "gwas.OA") := .(flip_strand(gwas.EA), flip_strand(gwas.OA))]
    map_ss[!(match), match := EA == gwas.EA | EA == gwas.OA]
    map_ss <- map_ss[(match)]
    if(nrow(map_ss) == 0) return(NULL)

    # Orient effect alleles
    map_ss[EA != gwas.EA, c("gwas.EA", "gwas.OA", "MAF", "effect") := .(gwas.OA, gwas.EA, 1 - MAF, -effect)]

    # Make sure ambiguous alleles are correctly oriented, if MAF has been provided
    map_ss[((ambig) | (gwas.ambig)) & ((EAF < 0.5 & MAF > 0.5) | (EAF > 0.5 & MAF < 0.5)),
            c("gwas.EA", "gwas.OA", "MAF", "effect") := .(gwas.OA, gwas.EA, 1 - MAF, -effect)]
    if(nrow(ind_ivs) == 0) return(NULL)

    # Orient to minor allele in INTERVAL
    map_ss[EAF > 0.5, c("EA", "OA", "effect", "EAF", "gwas.EA", "gwas.OA", "MAF", "gwas.effect") :=
                      .(OA, EA, -effect, 1-EAF, gwas.OA, gwas.EA, 1-MAF, -gwas.effect)]

    # Set up input data structures for co-localisation:
    pqtl_chr_stats <- list(
      pvalues = map_ss$P,
      beta = map_ss$effect,
      varbeta = (map_ss$se)^2,
      MAF = map_ss$EAF,
      N = 3301,
      type = "quant",
      sdY = 1,
      snp = map_ss[, paste(chr, pos, sep=":")]
    )

    # Input depends on whether the trait is case-control or continuous:
    if ("cases" %in% names(gwas_info)) {
      gwas_chr_stats <- list(
        pvalues = map_ss$gwas.P,
        beta = map_ss$gwas.effect,
        varbeta = (map_ss$gwas.se)^2,
        MAF = map_ss$MAF,
        N = gwas_info$samples,
        type = "cc",
        s = gwas_info$cases / gwas_info$samples,
        snp = map_ss[, paste(chr, pos, sep=":")]
      )
    } else {
      gwas_chr_stats <- list(
        pvalues = map_ss$gwas.P,
        beta = map_ss$gwas.effect,
        varbeta = (map_ss$gwas.se)^2,
        MAF = map_ss$MAF,
        N = gwas_info$samples,
        type = "quant",
        sdY = gwas_info$trait_sd,
        snp = map_ss[, paste(chr, pos, sep=":")]
      )
    }

    # Run coloc
    coloc <- coloc.abf(pqtl_chr_stats, gwas_chr_stats)

    # Collate posterior probabilities
    coloc_pp <- as.data.table(lapply(coloc$summary, `[`))

    # build return data.table
    test_info <- data.table(PRS = this_prs, GWAS = this_gwas, Aptamer = this_somaid, pQTL.chr = this_chr, pQTL.pos = this_pos)
    cbind(test_info, coloc_pp)
  }
}

# Call colocalisation
coloc[, colocalises := ((PP.H3.abf + PP.H4.abf) >= 0.9) & ((PP.H4.abf/PP.H3.abf) >= 3)]

# Add Gene
coloc[cis_info, on = .(Aptamer=SOMAMER_ID), Gene := Gene.Name]

fwrite(coloc, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/cis_iv_coloc.txt")

# Collate MR estimates
mr_collated <- copy(mr_summary)
mr_collated[prot_mr[Method == "(Intercept)"], on = .(PRS, GWAS, Gene), pleiotropy.pval := i.mr_pval]
mr_collated[coloc[,.(coloc=any(colocalises)), by=.(PRS, GWAS, Gene)], on = .(PRS, GWAS, Gene), coloc := i.coloc]
mr_collated <- merge(mr_collated, prot_mr, by = c("PRS", "GWAS", "Target", "UniProt", "Gene"), suffixes=c(".summary", ".protein"))
mr_collated <- merge(mr_collated, apt_mr, by = c("PRS", "GWAS", "Target", "UniProt", "Gene", "Method"))
mr_collated[cis_info, on = .(Aptamer=SOMAMER_ID), Aptamer := SeqId]
mr_collated[, Method := factor(Method, levels=c("IVW", "Median", "Weighted Median", "Weighted Mode", "MR Egger", "(Intercept)"))]
mr_collated <- mr_collated[order(Method)][order(Aptamer)][order(mr_pval.summary)]
mr_collated[, c("mr_se.summary", "mr_se.protein", "mr_se") := NULL]
mr_collated[, c("mr_estimate.summary", "mr_L95.summary", "mr_U95.summary", 
                "mr_estimate.protein", "mr_L95.protein", "mr_U95.protein", 
                "mr_estimate", "mr_L95", "mr_U95") :=
              .(exp(mr_estimate.summary), exp(mr_L95.summary), exp(mr_U95.summary), 
                exp(mr_estimate.protein), exp(mr_L95.protein), exp(mr_U95.protein),
                exp(mr_estimate), exp(mr_L95), exp(mr_U95))]

fwrite(mr_collated, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/collated_mr.txt")

# Collated instruments
ivs_collated <- copy(ind_ivs)
ivs_collated <- ivs_collated[coloc, on = .(PRS, GWAS, Gene, Aptamer, pQTL.chr, pQTL.pos)]
ivs_collated[cis_info, on = .(Aptamer=SOMAMER_ID), Aptamer := SeqId]
ivs_collated <- ivs_collated[order(pQTL.P)][order(Aptamer)]
ivs_collated <- ivs_collated[unique(mr_collated[,.(PRS, Target, UniProt, Gene, GWAS)]), on = .(PRS, Target, UniProt, Gene, GWAS)]
ivs_collated <- ivs_collated[, .(PRS, Target, UniProt, Gene, Gene.chr, Gene.start, 
  Aptamer, pQTL.chr, pQTL.pos, pQTL.EA, pQTL.OA, pQTL.EAF, pQTL.beta, pQTL.SE, pQTL.P,
  GWAS, gwas.EAF, gwas.beta, gwas.SE, gwas.P, nsnps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, colocalises)]

fwrite(ivs_collated, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/collated_ivs.txt")

# Plot IVs and mr estimates

# Average effects across aptamers (for IVs used in multiple aptamers)
gg_ivs <- ivs_collated[, .(pQTL.beta=mean(pQTL.beta), pQTL.SE=mean(pQTL.SE), pQTL.P=prot_pvalue(pQTL.P, pQTL.beta)),
                       by=.(PRS, GWAS, Gene, pQTL.chr, pQTL.pos, pQTL.EA, pQTL.OA, gwas.beta, gwas.SE, gwas.P)]
gg_ivs[, GWAS := fcase(GWAS == "T2DadjBMI", "T2D", GWAS == "StrokeIS", "IS", GWAS == "CAD", "CAD", GWAS == "CKD", "CKD")]

gg_mr <- unique(mr_collated[,.(PRS, GWAS, Gene, mr_estimate.summary, mr_L95.summary, mr_U95.summary)])
gg_mr <- gg_mr[unique(mr_collated[Method == "(Intercept)", .(PRS, GWAS, Gene, mr_estimate.protein, mr_L95.protein, mr_U95.protein)]), on = .(PRS, GWAS, Gene)]
setnames(gg_mr, c("PRS", "GWAS", "Gene", "OR", "L95", "U95", "Intercept", "Intercept.L95", "Intercept.U95"))
gg_mr[, GWAS := fcase(GWAS == "T2DadjBMI", "T2D", GWAS == "StrokeIS", "IS", GWAS == "CAD", "CAD", GWAS == "CKD", "CKD")]
gg_mr[, labtext := sprintf("%s ~ %s", Gene, GWAS)]
gg_mr[, labtext := factor(labtext, levels=unique(labtext))]
gg_ivs[gg_mr, on = .(PRS, GWAS, Gene), labtext := labtext]

# Determine the plot limits
plotlim <- gg_ivs[, .(xmin = min(pQTL.beta - pQTL.SE), xmax = max(pQTL.beta + pQTL.SE),
                      ymin = min(gwas.beta - gwas.SE), ymax = max(gwas.beta + gwas.SE))]

# Make sure to include intercepts
plotlim[xmax < 0, xmax := 0]
plotlim[xmin > 0, xmin := 0]
plotlim[ymax < 0, ymax := 0]
plotlim[ymin > 0, ymin := 0]

# Expand by 5% as per ggplot defaults
expand <- plotlim[, .(xexpand = abs((xmax - xmin))*0.05, yexpand = abs((ymax - ymin))*0.05)]
plotlim[, c("xmin", "xmax", "ymin", "ymax") := .(xmin - expand$xexpand, xmax + expand$xexpand, ymin - expand$yexpand, ymax + expand$yexpand)]

# Now build polygons for 95% CIs. This table is a series of x and y-coordinates for drawing
# the polygon for each protein to disease pair. We need to work out where the polygon intersects
# the edges of the plot so we know where to define the polygon edges.
gg_dr_ci95 <- foreach(testIdx = gg_mr[,.I], .combine=rbind) %do% {
  mrtest <- gg_mr[testIdx]
  intercept.l95 <- mrtest[, log(Intercept.L95)]
  intercept.u95 <- mrtest[, log(Intercept.U95)]
  slope.l95 <- mrtest[, log(L95)]
  slope.u95 <- mrtest[, log(U95)]

  xlim <- plotlim[, c(xmin, xmax)]
  ylim <- plotlim[, c(ymin, ymax)]

  if (intercept.l95 < ylim[1]) intercept.l95 <- ylim[1]
  if (intercept.u95 > ylim[2]) intercept.u95 <- ylim[2]

  poly <- line_ci95_poly(intercept.l95, intercept.u95, slope.l95, slope.u95, xlim, ylim)
  cbind(gg_mr[testIdx,.(labtext)], poly)
}

# plot
g <- ggplot(gg_ivs) +
  aes(x=pQTL.beta, xmin=pQTL.beta - pQTL.SE, xmax=pQTL.beta + pQTL.SE,
      y=exp(gwas.beta), ymin=exp(gwas.beta - gwas.SE), ymax=exp(gwas.beta + gwas.SE)) +
  geom_polygon(data=gg_dr_ci95, inherit.aes=FALSE, aes(x=x, y=exp(y)), fill="#fed976", color=NA, size=0) +
  geom_hline(yintercept=1, linetype=3, color="#bdbdbd", size=0.2) +
  geom_vline(xintercept=0, linetype=3, color="#bdbdbd", size=0.2) +
  geom_abline(data=gg_mr, aes(intercept=Intercept, slope=-(1 - OR)), color="#e31a1c", linetype=2, size=0.3) +
  geom_errorbarh(height=0, size=0.2) + geom_errorbar(width=0, size=0.2) +
  geom_point(shape=23, size=0.6, fill="white", stroke=0.4) +
  facet_wrap(~ labtext, nrow=3, dir = "h") +
  scale_x_continuous(name="Effect of pQTL on protein", expand=c(0,0)) +
  scale_y_continuous(name="OR for pQTL on disease", expand=c(0,0)) +
  theme_bw() +
  theme(axis.title=element_text(size=7), axis.text=element_text(size=5), panel.grid=element_blank(),
        strip.text=element_blank(), strip.background=element_blank(), panel.spacing=unit(0.5, "mm"),
        axis.ticks=element_line(size=0.3), plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"))
ggsave(g, width=109, height=46, units="mm", file="analyses/pub/cardiometabolic_proteins/review4/mr_dose_response.pdf")

