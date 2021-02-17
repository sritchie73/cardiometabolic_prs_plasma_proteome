library(data.table)
library(foreach)
library(doMC)
library(coloc)
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

# Other options
biallelic_only <- TRUE
snps_only <- TRUE

# Load cis IVs
cis_ivs <- fread("analyses/pub/cardiometabolic_proteins/review2/cis_ivs.txt")

# Only test colocalisation if GWAS p-value < 1-e6
cis_ivs <- cis_ivs[gwas.P < 1e-6]

# Load protein information to get SOMAMER_ID
apt_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
apt_info <- apt_info[,.(Aptamer=SeqId, SOMAMER_ID)]

# Test colocalisation for each of these instruments
gwass <- unique(cis_ivs[,.(PRS, GWAS)])
gwass[, gwas_dir := sprintf("analyses/mendelian_randomisation/%s/gwas_summary_stats/%s", PRS, GWAS)]
coloc <- foreach(gwas_idx = gwass[,.I], .combine=rbind) %do% {
  this_prs <- gwass[gwas_idx, PRS]
  this_gwas <- gwass[gwas_idx, GWAS]

  # Determine windows we need to load for this GWAS (200Kb around each pQTL)
  windows <- cis_ivs[PRS == this_prs & GWAS == this_gwas, .(
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
  this_cis_ivs <- cis_ivs[PRS == this_prs & GWAS == this_gwas]
  foreach(test_idx = this_cis_ivs[,.I], .combine=rbind) %dopar% {
    this_chr <- this_cis_ivs[test_idx, pQTL.chr]
    this_pos <- this_cis_ivs[test_idx, pQTL.pos]
    window_start <- this_pos - 200000
    window_end <- this_pos + 200000
    this_aptid <- this_cis_ivs[test_idx, Aptamer]
    this_somaid <- apt_info[Aptamer == this_aptid, SOMAMER_ID]

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
    if(nrow(cis_ivs) == 0) return(NULL)

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
    test_info <- data.table(PRS = this_prs, GWAS = this_gwas, Aptamer = this_aptid, pQTL.chr = this_chr, pQTL.pos = this_pos)
    cbind(test_info, coloc_pp)
  }
}

# Call colocalisation
coloc[, colocalises := ((PP.H3.abf + PP.H4.abf) >= 0.9) & ((PP.H4.abf/PP.H3.abf) >= 3)]

fwrite(coloc, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/cis_iv_coloc.txt")





