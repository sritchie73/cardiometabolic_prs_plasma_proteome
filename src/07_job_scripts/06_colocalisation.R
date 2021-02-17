library(data.table)
library(foreach)
library(coloc)
library(doMC)
source("src/utilities/flip_strand.R")

# Determine GRS we're working with
if (!exists("grs")) {
  view_file <- commandArgs(trailingOnly=TRUE)[1]
  view <- fread(view_file)
  task <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  grs <- view[task, GRS_name]
}

# Determine number of CPUs available
ncores <- as.numeric(system("echo $SLURM_CPUS_ON_NODE", intern=TRUE))
ncores[is.na(ncores)] <- 1
registerDoMC(1)

# Paths:
out_dir <- sprintf("analyses/mendelian_randomisation/%s", grs)
pqtl_dir <- "data/full_pQTL_summary_stats"
gwas_dir <- sprintf("analyses/mendelian_randomisation/%s/gwas_summary_stats", grs)
snp_dir <- "data/INTERVAL/reference_files/imputed_genotypes"

# Check if there is anything to run:
if (!file.exists(sprintf("%s/instruments.tsv", out_dir))) {
  quit(save="no", status=0)
}
if (dir.exists(sprintf("%s/coloc", out_dir))) {
  quit(save="no", status=0)
}

# load instruments table
instruments <- fread(sprintf("%s/instruments.tsv", out_dir))

# Order by chromosome to reduce number of reads of snpstats
instruments <- instruments[order(chr)]

# Run colocalisation for each instrument:
last_chr <- -1
snpstats <- data.table()
coloc_dt <- foreach(rIdx = 1:nrow(instruments), .combine=rbind, .errorhandling="remove") %dopar% {
  aptamer <- instruments[rIdx, SOMAMER_ID]
  gwas <- instruments[rIdx, GWAS]
  instrument <- instruments[rIdx, ALT_ID]
  pqtl_chr <- instruments[rIdx, chr]
  pqtl_pos <- instruments[rIdx, pos]
  gwas_info <- fread(sprintf("%s/%s/info.txt", gwas_dir, gwas))

  # Only test colocalisation if GWAS p-value < 1-e6
  if (instruments[rIdx, P.gwas >= 1e-6]) {
    coloc_pp <- data.table(nsnps=NA_integer_, PP.H0.abf=NA_real_, PP.H1.abf=NA_real_, 
                           PP.H2.abf=NA_real_, PP.H3.abf=NA_real_, PP.H4.abf=NA_real_)
    coloc_pp <- cbind(instruments[rIdx, .(GWAS, SOMAMER_ID, ALT_ID)], coloc_pp)
    return(coloc_pp)
  }

  # Determine window:
  dist <- 200000
  window_start <- pmax(pqtl_pos - dist, 0)
  window_end <- pqtl_pos + dist

  # Load snpstats so we can get MAF - preserve between iterations if on the 
  # same chromosome.
  if (pqtl_chr != last_chr) {
    last_chr <- pqtl_chr
    snpstats <<- fread(sprintf("%s/impute_%s_interval.snpstats", snp_dir, pqtl_chr)) 
  }

  # Load pQTL summary statistics for this aptamer and filter variants 200kb up or downstream
  pqtl_chr_stats <- fread(sprintf("%s/%s/%s_chrom_%s_meta_1.tbl.gz", pqtl_dir, aptamer, aptamer, pqtl_chr))
	pqtl_chr_stats <- pqtl_chr_stats[, .(pos=position, EA=toupper(Allele1), OA=toupper(Allele2),
																				 effect=Effect, se=StdErr, P=10^`log(P)`)]
  pqtl_chr_stats <- pqtl_chr_stats[(pos >= (window_start)) & (pos <= (window_end))]

  # Add MAF of effect allele
  pqtl_chr_stats[snpstats, on = .(pos=position, EA=minor_allele, OA=major_allele), c("EAF", "RSID") := .(i.MAF, RSID)]
  pqtl_chr_stats[snpstats, on = .(pos=position, EA=major_allele, OA=minor_allele), c("EAF", "RSID") := .(1 - i.MAF, RSID)]

	# Filter out variants that can't be reliably mapped due to strand ambiguity
	pqtl_chr_stats <- pqtl_chr_stats[EA != flip_strand(OA)]

  # Load GWAS summary stats and filter to variants 200kb up or downstream of the pQTL.
	gwas_chr_stats <- fread(sprintf("%s/%s/chr%s.txt.gz", gwas_dir, gwas, pqtl_chr))
  gwas_chr_stats <- gwas_chr_stats[, .(pos, EA, OA, EAF=MAF, effect, se, P)]
  gwas_chr_stats <- gwas_chr_stats[(pos >= (pqtl_pos - 200000)) & (pos <= (pqtl_pos + 200000))]

	# Filter out variants that can't be reliably mapped due to strand ambiguity
	gwas_chr_stats <- gwas_chr_stats[EA != flip_strand(OA)]

	# Get pQTLs and proxies that can be directly matched between datasets:
	strand_ok_col_ok <- merge(pqtl_chr_stats, gwas_chr_stats,
														by=c("pos", "EA", "OA"), suffixes=c(".pqtl", ".gwas"))

	# Handle cases where the effect allele is the other allele in the GWAS dataset
	strand_ok_col_swap <- merge(pqtl_chr_stats, gwas_chr_stats,
															by.x=c("pos", "EA", "OA"), by.y=c("pos", "OA", "EA"),
															suffixes=c(".pqtl", ".gwas"))
	strand_ok_col_swap[, c("effect.gwas", "EAF.gwas") := .(-effect.gwas, 1 - EAF.gwas)]

	# Now get ones that can be directly matched after flipping the strand:
	pqtl_chr_stats[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]
	strand_flip_col_ok <- merge(pqtl_chr_stats, gwas_chr_stats,
															by=c("pos", "EA", "OA"), suffixes=c(".pqtl", ".gwas"))
	strand_flip_col_ok[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

	# Strand flip, effect allele = other allele
	strand_flip_col_swap <- merge(pqtl_chr_stats, gwas_chr_stats,
																by.x=c("pos", "EA", "OA"), by.y=c("pos", "OA", "EA"),
																suffixes=c(".pqtl", ".gwas"))
	strand_flip_col_swap[, c("effect.gwas", "EAF.gwas") := .(-effect.gwas, 1 - EAF.gwas)]
	strand_flip_col_swap[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

	# Reset pqtl_chr_stats to original strand
	pqtl_chr_stats[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

	# Alleles always oriented to effect allele and strand of pQTL study:
	combined <- rbind(strand_ok_col_ok, strand_ok_col_swap,
										strand_flip_col_ok, strand_flip_col_swap)

  # If MAF info is missing from the GWAS summary stats, use the MAF from INTERVAL
  combined[is.na(EAF.gwas), EAF.gwas := EAF.pqtl]

  # Orient to minor allele in INTERVAL
  combined[EAF.pqtl > 0.5, 
    c("EAF.gwas", "EAF.pqtl", "EA", "OA", "effect.gwas", "effect.pqtl") :=
    .(1 - EAF.gwas, 1 - EAF.pqtl, OA, EA, -effect.gwas, -effect.pqtl)]

  # Fix identifiers
  combined[is.na(RSID), RSID := "."]
  combined[, ALT_ID := paste0("chr", pqtl_chr, ":", pos, ":", EA, ":", OA)]

  # Set up input data structures for co-localisation:
  pqtl_chr_stats <- list(
    pvalues = combined$P.pqtl,
    beta = combined$effect.pqtl,
    varbeta = (combined$se.pqtl)^2,
    MAF = combined$EAF.pqtl,
    N = 3301,
    type = "quant",
    sdY = 1,
    snp = combined$ALT_ID
  )

  # Input depends on whether the trait is case-control or continuous:
  if ("cases" %in% names(gwas_info)) {
		gwas_chr_stats <- list(
			pvalues = combined$P.gwas,
			beta = combined$effect.gwas,
			varbeta = (combined$se.gwas)^2,
			MAF = combined$EAF.gwas,
			N = gwas_info$samples,
			type = "cc",
			s = gwas_info$cases / gwas_info$samples,
			snp = combined$ALT_ID
		)
  } else {
		gwas_chr_stats <- list(
			pvalues = combined$P.gwas,
			beta = combined$effect.gwas,
			varbeta = (combined$se.gwas)^2,
			MAF = combined$EAF.gwas,
			N = gwas_info$samples,
			type = "quant",
			sdY = gwas_info$trait_sd,
			snp = combined$ALT_ID
		)
  }
 
  # Run co-localisation
  coloc <- coloc.abf(pqtl_chr_stats, gwas_chr_stats)

  # Output details:
  coloc_dir <- sprintf("%s/coloc/%s/%s", out_dir, aptamer, gwas)
  dir.create(coloc_dir, recursive=TRUE, showWarnings=FALSE)
 
  details <- as.data.table(coloc$results)
  fwrite(details, sep="\t", quote=FALSE, 
    file=sprintf("%s/%s_coloc_results.tsv", coloc_dir, pqtl_pos))
  fwrite(combined, sep="\t", quote=FALSE, 
    file=sprintf("%s/%s_coloc_input.tsv", coloc_dir, pqtl_pos))

  # Collate posterior probabilities
  coloc_pp <- as.data.table(lapply(coloc$summary, `[`))
  coloc_pp <- cbind(instruments[rIdx, .(GWAS, SOMAMER_ID, ALT_ID)], coloc_pp)
  return(coloc_pp)
}

coloc_dt[, colocalises := ((PP.H3.abf + PP.H4.abf) >= 0.9) & ((PP.H4.abf/PP.H3.abf) >= 3)]
coloc_dt[is.na(colocalises), colocalises := FALSE]

instruments[coloc_dt, on=.(GWAS, SOMAMER_ID, ALT_ID),
  c("N_coloc_snps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "colocalises") :=
  .(i.nsnps, i.PP.H0.abf, i.PP.H1.abf, i.PP.H2.abf, i.PP.H3.abf, i.PP.H4.abf, i.colocalises)]
 
fwrite(instruments, sep="\t", quote=FALSE, file=sprintf("%s/instruments.tsv", out_dir))

# Create LoucsZoom style plots afterwards. To do this, we need to know the LD
# with the instrument in each 1MB window. Instead of calculating LD in INTERVAL
# (slow) we'll query the NIH NCI LDlink API. 
token <- "eb688cc2c097"
to_plot <- instruments[!is.na(PP.H4)]
foreach(iv = unique(to_plot$ALT_ID), .errorhandling="pass") %dopar% {
  iv_pos <- to_plot[ALT_ID == iv, pos][1]
  iv_chr <- to_plot[ALT_ID == iv, chr][1]
  coord_id <- paste0("chr", iv_chr, ":", iv_pos)
  window_start <- pmax(iv_pos - 200000, 0)
  window_end <- iv_pos + 200000

  # Pull down information about proxies for this variant. Note, LD proxy fortunately
  # uses hg19.
  tryCatch({
		uri <- sprintf("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=%s&pop=GBR&r2_d=r2&token=%s", coord_id, token)
		cmd <- sprintf("curl -k -X GET '%s'", uri)
		ld <- fread(cmd=cmd)
		ld <- ld[, .(pos=as.numeric(gsub("chr[0-9]*:", "", Coord)), r2=R2)]
  }, error=function(e) {
    # If we reach this block it means the LDlink request has failed.
    # This should only happen if:
    #  (1) The variant is multiallelic
    #  (2) The variant is not a SNP or INDEL
    cat(sprintf("LDlink lookup failed for variant %s, calculating LD in INTERVAL instead\n", iv), file=stderr())

    # Directory to put temporary files in
    ld_dir <- sprintf("%s/coloc", out_dir)

		# Get r2 with lead pQTL for locus-zoom style plots
		# First use bgenix to extract 200kb window
		cmd <- "bgenix"
		cmd <- paste(cmd, sprintf("-g data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/impute_%s_interval.bgen", iv_chr))
		if (iv_chr < 10) { # %01d doesnt seem to work here??
			cmd <- paste(cmd, sprintf("-incl-range 0%s:%s-%s", iv_chr, window_start, window_end))
		} else {
			cmd <- paste(cmd, sprintf("-incl-range %s:%s-%s", iv_chr, window_start, window_end))
		}
		cmd <- paste(cmd, sprintf("> %s/%s.window.bgen", ld_dir, iv_pos))
		system(command=cmd, wait=TRUE)

		# Then extract pQTL into separate BGEN file so we can compute its LD with other
		# variants in the window using qctool2 (ldstore is too slow because it calculates
		# all pairwise in window).
		pQTL_var <- snpstats[position == iv_pos, .(SNPID, RSID, chromosome, position, A_allele, B_allele)]
		pQTL_var[, chromosome := ifelse(iv_chr < 10, paste0("0", iv_chr), iv_chr)]
		fwrite(pQTL_var, sep=" ", quote=FALSE,
					 file=sprintf("%s/%s.var.txt", ld_dir, iv_pos))

		cmd <- "qctool2"
		cmd <- paste(cmd, sprintf("-g %s/%s.window.bgen", ld_dir, iv_pos))
		cmd <- paste(cmd, sprintf("-og %s/%s.bgen", ld_dir, iv_pos))
		cmd <- paste(cmd, sprintf("-incl-variants %s/%s.var.txt", ld_dir, iv_pos))
		system(cmd, wait=TRUE)

		# Extract variant set to calculate LD for
		other_var <- snpstats[position %in% combined$pos][position != iv_pos]
		other_var <- other_var[, .(SNPID, RSID, chromosome, position, A_allele, B_allele)]
		other_var[, chromosome := ifelse(iv_chr < 10, paste0("0", iv_chr), as.character(iv_chr))]
		fwrite(other_var, sep=" ", quote=FALSE,
					 file=sprintf("%s/%s.var2.txt", ld_dir, iv_pos))

		cmd <- "qctool2"
		cmd <- paste(cmd, sprintf("-g %s/%s.window.bgen", ld_dir, iv_pos))
		cmd <- paste(cmd, sprintf("-og %s/%s.other.bgen", ld_dir, iv_pos))
		cmd <- paste(cmd, sprintf("-incl-variants %s/%s.var2.txt", ld_dir, iv_pos))
		system(cmd, wait=TRUE)

		# Use qctool2 to compute LD between pQTL and all other variants in window:
		cmd <- "qctool2"
		cmd <- paste(cmd, sprintf("-g %s/%s.bgen", ld_dir, iv_pos))
		cmd <- paste(cmd, "-s data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/interval.samples")
		cmd <- paste(cmd, sprintf("-compute-ld-with %s/%s.window.bgen", ld_dir, iv_pos))
		cmd <- paste(cmd, "data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/interval.samples")
		cmd <- paste(cmd, sprintf("-old %s/%s.ld.sqlite", ld_dir, iv_pos))
		cmd <- paste(cmd, "-min-r2 0.2")
		system(cmd, wait=TRUE)

		# Load LD
		cmd <- "sqlite3 -header -separator $'\t'"
		cmd <- paste(cmd, sprintf("%s/%s.ld.sqlite", ld_dir, iv_pos))
		cmd <- paste(cmd, '"SELECT * FROM ld_analysis1View"')
		ld <- fread(cmd=cmd)
    ld <- ld[, .(pos=variant1_position, r2=dosage_r2)]

		# Clean up
		system(sprintf("rm %s/*.*", ld_dir), wait=TRUE)
  })
   
  foreach(aptamer = unique(to_plot[ALT_ID == iv, SOMAMER_ID]), .errorhandling="pass") %do% {
    foreach(gwas = to_plot[ALT_ID == iv & SOMAMER_ID == aptamer, GWAS], .errorhandling="pass") %do% {
      # load data
      coloc_dir <- sprintf("%s/coloc/%s/%s", out_dir, aptamer, gwas)
      coloc_input <- fread(sprintf("%s/%s_coloc_input.tsv", coloc_dir, iv_pos))

      # Build table of LD buckets
      ld <- rbind(ld, coloc_input[!ld, on=.(pos), .(pos, r2=0)])
      ld[r2 >= 0.8 & r2 < 1, r2 := 0.8]
      ld[r2 >= 0.6 & r2 < 0.8, r2 := 0.6]
      ld[r2 >= 0.4 & r2 < 0.6, r2 := 0.4]
      ld[r2 >= 0.2 & r2 < 0.4, r2 := 0.2]
      ld[r2 < 0.2, r2 := 0]
      ld <- ld[pos != iv_pos]

      # Generate Locus-Zoom style plots:
      pdf(sprintf("%s/%s_coloc.pdf", coloc_dir, iv_pos), width=3.6, height=7.2)
      par(mfrow=c(2,1))
      xlim <- c(window_start, window_end) + c(-400000 * 0.01, 400000 * 0.01)
      ylim <- c(0, max(-log10(coloc_input$P.pqtl)))
      ylim <- ylim + c(0, ylim[2] * 0.01)
      plot(0, type="n", xlim=xlim, ylim=ylim, 
           xlab=paste("Position on chromosome", iv_chr),
           ylab="-log10 P-value in pQTL study")
      coloc_input[ld[r2 == 0.0], on = .(pos), nomatch=0, points(pos, -log10(P.pqtl), pch=21, bg="#000080", col="black")]
      coloc_input[ld[r2 == 0.2], on = .(pos), nomatch=0, points(pos, -log10(P.pqtl), pch=21, bg="#87cefa", col="black")]
      coloc_input[ld[r2 == 0.4], on = .(pos), nomatch=0, points(pos, -log10(P.pqtl), pch=21, bg="#00ff00", col="black")]
      coloc_input[ld[r2 == 0.6], on = .(pos), nomatch=0, points(pos, -log10(P.pqtl), pch=21, bg="#ffa500", col="black")]
      coloc_input[ld[r2 == 0.8], on = .(pos), nomatch=0, points(pos, -log10(P.pqtl), pch=21, bg="#ff0000", col="black")]
      coloc_input[ld[r2 == 1], on = .(pos), nomatch=0, points(pos, -log10(P.pqtl), pch=21, bg="#7d26cd", col="black")]
      coloc_input[pos == iv_pos, points(pos, -log10(P.pqtl), pch=23, bg="#7d26cd", col="black")]

      ylim <- c(0, max(-log10(coloc_input$P.gwas)))
      ylim <- ylim + c(0, ylim[2] * 0.01)
      plot(0, type="n", xlim=xlim, ylim=ylim, 
           xlab=paste("Position on chromosome", iv_chr),
           ylab="-log10 P-value in GWAS")
      coloc_input[ld[r2 == 0.0], on = .(pos), nomatch=0, points(pos, -log10(P.gwas), pch=21, bg="#000080", col="black")]
      coloc_input[ld[r2 == 0.2], on = .(pos), nomatch=0, points(pos, -log10(P.gwas), pch=21, bg="#87cefa", col="black")]
      coloc_input[ld[r2 == 0.4], on = .(pos), nomatch=0, points(pos, -log10(P.gwas), pch=21, bg="#00ff00", col="black")]
      coloc_input[ld[r2 == 0.6], on = .(pos), nomatch=0, points(pos, -log10(P.gwas), pch=21, bg="#ffa500", col="black")]
      coloc_input[ld[r2 == 0.8], on = .(pos), nomatch=0, points(pos, -log10(P.gwas), pch=21, bg="#ff0000", col="black")]
      coloc_input[ld[r2 == 1], on = .(pos), nomatch=0, points(pos, -log10(P.gwas), pch=21, bg="#7d26cd", col="black")]
      coloc_input[pos == iv_pos, points(pos, -log10(P.gwas), pch=23, bg="#7d26cd", col="black")]
      dev.off()
    }
  }
}

