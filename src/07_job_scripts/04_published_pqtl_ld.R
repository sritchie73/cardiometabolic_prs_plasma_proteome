library(data.table)
library(foreach)
library(doMC)

args <- commandArgs(trailingOnly=TRUE)
nCores <- args[1]
registerDoMC(nCores)

out_dir <- "analyses/mendelian_randomisation/pqtls/conditional_pQTLs_tag_LD"

# Get list of conditionally independent pQTLs
if (file.exists("analyses/mendelian_randomisation/pqtls/conditionally_independent_pQTLs.txt")) {
  cond_pQTLs <- fread("analyses/mendelian_randomisation/pqtls/conditionally_independent_pQTLs.txt")
} else {
	cond_pQTLs <- fread("data/raw/INTERVAL/pQTLs.txt")
	cond_pQTLs <- cond_pQTLs[, .(SOMAMER_ID=`SOMAmer ID`, rsID=`Conditional variant`, sentinel=`Sentinel variant`)]
	pQTL_info <- fread("data/raw/INTERVAL/pQTL_sentinel_info.txt")
	cond_pQTLs[pQTL_info, on=.(SOMAMER_ID=`SOMAmer ID`, sentinel=`Sentinel variant*`), type := `cis/ trans`]

	# Get chromosomal locations for each conditional pQTL (data in loaded file not reliable,
	# same rsID is listed as having multiple genomic positions, likely corresponds to sentinel variant).
	cond_pQTLs[grepl(":", rsID), position := as.numeric(gsub(".*:", "", rsID))]
	cond_pQTLs[grepl(":", rsID), chromosome := as.numeric(gsub("chr", "", gsub(":.*", "", rsID)))]
	cond_pQTLs[grepl(":", rsID), rsID := "."]

	for (chrIdx in 1:22) {
		chr_stats <- fread(sprintf("data/raw/INTERVAL/reference_files/impute_%s_interval.snpstats", chrIdx),
											 colClasses=c("chromosome"="numeric", "position"="numeric"))
		cond_pQTLs[chr_stats[RSID != "."], on=.(rsID=RSID), 
							 c("chromosome", "position", "A_allele", "B_allele") :=
							 .(i.chromosome, i.position, i.A_allele, i.B_allele)]
		cond_pQTLs[chr_stats[RSID == "."], on=.(chromosome, position),
							 c("rsID", "A_allele", "B_allele") :=
							 .(RSID, i.A_allele, i.B_allele)]
		dir.create(sprintf("%s/chr%s/", out_dir, chrIdx), showWarnings=FALSE)
	}
	fwrite(cond_pQTLs, sep="\t", quote=FALSE,
				 file="analyses/mendelian_randomisation/pqtls/conditionally_independent_pQTLs.txt")
}

# Filter to unique variants - we will run all cis and trans pQTLs, because the
# cis-pQTL run only looked at SNPs with P < 1e-4 with unambiguous alleles, where
# here we need to make sure the conditional variants have their LD calculated,
# so we can identify proxies in the case they had ambiguous alleles.
cond_pQTLs <- unique(cond_pQTLs[, .(rsID, chromosome, position, A_allele, B_allele)])

# For each pQTL, get all variants with r2 > 0.8 within a 250KB window
foreach (rIdx = 1:nrow(cond_pQTLs), .errorhandling="remove") %dopar% {
  chr <- as.numeric(cond_pQTLs[rIdx, chromosome])
  pos <- cond_pQTLs[rIdx, position]
  window_start <- ifelse(pos < 250000, 0, pos - 250000)
  window_end <- pos + 250000

  # First use bgenix to extract 250kb window
  cmd <- "bgenix"
  cmd <- paste(cmd, sprintf("-g data//INTERVAL/impute_%s_interval.bgen", chr))
  if (chr < 10) { # %01d doesnt seem to work here??
    cmd <- paste(cmd, sprintf("-incl-range 0%s:%s-%s", chr, window_start, window_end))
  } else {
    cmd <- paste(cmd, sprintf("-incl-range %s:%s-%s", chr, window_start, window_end))
  }
  cmd <- paste(cmd, sprintf("> %s/chr%s/%s.250kb_window.bgen", out_dir, chr, pos))
  system(command=cmd, wait=TRUE)

  # Then extract pQTL into separate BGEN file so we can compute its LD with other
  # variants in the window using qctool2 (ldstore is too slow because it calculates
  # all pairwise in window).
  pQTL_var <- cond_pQTLs[rIdx, .(SNPID=rsID, RSID=rsID, chromosome, position, A_allele, B_allele)]
  pQTL_var[, chromosome := ifelse(chr < 10, paste0("0", chr), chr)]
  fwrite(pQTL_var, sep=" ", quote=FALSE, 
         file=sprintf("%s/chr%s/%s.var.txt", out_dir, chr, pos))

  cmd <- "qctool2"
  cmd <- paste(cmd, sprintf("-g %s/chr%s/%s.250kb_window.bgen", out_dir, chr, pos))
  cmd <- paste(cmd, sprintf("-og %s/chr%s/%s.bgen", out_dir, chr, pos))
  cmd <- paste(cmd, sprintf("-incl-variants %s/chr%s/%s.var.txt", out_dir, chr, pos))
  system(cmd, wait=TRUE)

  # Use qctool2 to compute LD between pQTL and all other variants in 250KB window:
  cmd <- "qctool2"
  cmd <- paste(cmd, sprintf("-g %s/chr%s/%s.bgen", out_dir, chr, pos))
  cmd <- paste(cmd, "-s data//INTERVAL/interval.samples")
  cmd <- paste(cmd, sprintf("-compute-ld-with %s/chr%s/%s.250kb_window.bgen", out_dir, chr, pos))
  cmd <- paste(cmd, "data//INTERVAL/interval.samples")
  cmd <- paste(cmd, sprintf("-old %s/chr%s/%s.ld.sqlite", out_dir, chr, pos))
  cmd <- paste(cmd, "-min-r2 0.8")
  system(cmd, wait=TRUE)

  # Extract results into a table mimicing the output of the cis-pQTL LD tables
  cmd <- "sqlite3 -header  -separator $'\t'"
  cmd <- paste(cmd, sprintf("%s/chr%s/%s.ld.sqlite", out_dir, chr, pos))
  cmd <- paste(cmd, '"SELECT * FROM ld_analysis1View"')
  ld <- fread(cmd=cmd)
  ld <- ld[, .(chromosome=variant1_chromosome,
               RSID1=variant2_rsid, position1=variant2_position,
               RSID2=variant1_rsid, position2=variant1_position,
               correlation=dosage_r, r2=dosage_r2)]
  ld <- ld[RSID1 != RSID2]
  fwrite(ld, sep="\t", quote=FALSE, file=sprintf("%s/chr%s/%s.ld", out_dir, chr, pos))

  # Clean up
  system(sprintf("rm %s/chr%s/%s.250kb_window.bgen", out_dir, chr, pos))
  system(sprintf("rm %s/chr%s/%s.var.txt", out_dir, chr, pos))
  system(sprintf("rm %s/chr%s/%s.bgen", out_dir, chr, pos))
  system(sprintf("rm %s/chr%s/%s.ld.sqlite", out_dir, chr, pos))
}

