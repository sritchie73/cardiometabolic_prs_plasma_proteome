library(data.table)
library(foreach)
library(doMC)
source("src/utilities/flip_strand.R")

nCores <- commandArgs(trailingOnly=TRUE)[1]
registerDoMC(nCores)

out_dir <- "analyses/mendelian_randomisation/pqtls/SOMAMERs_cis_LD"

# Get TSS info about all proteins with cis-pQTLs
info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
cis <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")

# Get relevant columns for filtering variants from BGEN files
to_filter <- foreach(chrIdx = 1:22, .combine=rbind) %dopar% {
	chr_stats <- fread(sprintf("data/INTERVAL/reference_files/imputed_genotypes/impute_%s_interval.snpstats", chrIdx))
	chr_stats <- rbind(
		chr_stats[cis, on = .(chromosome=chr, position=pos, A_allele=EA, B_allele=OA), nomatch=0],
		chr_stats[cis, on = .(chromosome=chr, position=pos, A_allele=OA, B_allele=EA), nomatch=0])
	chr_stats <- chr_stats[, .(SOMAMER_ID, RSID, position, chromosome, A_allele, B_allele)]
	return(chr_stats)
}

# Ignore ambiguous variants, these may confuse the calculations
# to_filter <- to_filter[A_allele != flip_strand(B_allele)] 

# Filter to aptamers with > 1 pqtl
n_pQTLs <- to_filter[,.N,by=.(SOMAMER_ID, chromosome)]

# Split out protein information to gene level information
info <- info[Type == "Protein", .(SOMAMER_ID, Gene.Name, chr, start)]
info <- info[chr != "" & start != ""]
locmap <- info[,.(
  Gene=strsplit(Gene.Name, "\\|")[[1]],
  chr=strsplit(chr, "\\|")[[1]],
  start=strsplit(start, "\\|")[[1]]),
  by=.(SOMAMER_ID)]
locmap <- locmap[, .(start = strsplit(start, "\\;")[[1]]), by = .(SOMAMER_ID, Gene, chr)] # IGK has multiple start positions
locmap <- locmap[!(chr %in% c("X", "Y"))] # no genotype data for X and Y chr, thus no cis-pQTLs
locmap[, start := as.integer(start)]
locmap[, chr := as.integer(chr)]

# Filter to genes with more than 1 cis-PQTL
locmap <- locmap[n_pQTLs[N > 1], on = .(SOMAMER_ID, chr=chromosome)]

# Calculate window to for which to consider for LD
locmap <- locmap[, .(window_start = min(start) - 1e6, 
                     window_end = max(start) + 1e6),
                 by=.(SOMAMER_ID, chr, N)]
locmap[window_start < 0, window_start := 0]

# Distribute long-running LD jobs (many pQTLs) evenly across cores.
# %dopar% assigns each element a new process until it runs out of 
# cores, then starts again (i.e. for 6 cores, elements 1, 7, 13, etc.
# will all be on core 1, elements 2, 8, 14 etc will be on core 2,
# and so on). So we can load balance by sorting by number of pQTLs
# to compute
locmap <- locmap[order(N)]

# Calculate LD for cis-pQTLs for each protein with > 1 pQTL
foreach (rIdx = 1:nrow(locmap), .errorhandling="remove") %dopar% {
  aptamer <- locmap[rIdx, SOMAMER_ID]
  chr <- as.numeric(locmap[rIdx, chr])
  window_start <- locmap[rIdx, window_start]
  window_end <- locmap[rIdx, window_end]

  # Write out variant file to use for filtering in qctool2
  prot_vars <- to_filter[SOMAMER_ID == aptamer] 
  prot_vars <- prot_vars[,.(SNPID=RSID, RSID=RSID, chromosome, position, A_allele, B_allele)]
  if (chr < 10) { # sprintf %01d does not work here?
    prot_vars[, chromosome := paste0("0", chromosome)]
  }
  fwrite(prot_vars, sep=" ", quote=FALSE, file=sprintf("%s/%s_%s.vars.txt", out_dir, aptamer, chr))
 
  # First use bgenix to extract 1MB window, to reduce memory usage and  
  # speed up variant extract by qctool2
  cmd <- "bgenix"
  cmd <- paste(cmd, sprintf("-g data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/impute_%s_interval.bgen", chr))
  if (chr < 10) { # %01d doesnt seem to work here??
		cmd <- paste(cmd, sprintf("-incl-range 0%s:%s-%s", chr, window_start, window_end))
  } else {
		cmd <- paste(cmd, sprintf("-incl-range %s:%s-%s", chr, window_start, window_end))
  }
  cmd <- paste(cmd, sprintf("> %s/%s_%s.cis_window.bgen", out_dir, aptamer, chr))
  system(command=cmd, wait=TRUE)

  # Run qctool2 command extracting the cis-pQTL variants into new BGEN file
  cmd <- "qctool2"
  cmd <- paste(cmd, sprintf("-g %s/%s_%s.cis_window.bgen", out_dir, aptamer, chr))
  cmd <- paste(cmd, sprintf("-og %s/%s_%s.bgen", out_dir, aptamer, chr))
  cmd <- paste(cmd, sprintf("-incl-variants %s/%s_%s.vars.txt", out_dir, aptamer, chr))
  system(cmd, wait=TRUE)

  # Then use LDstore to calculate LD between variants
  cmd <- "ldstore"
  cmd <- paste(cmd, sprintf("--bcor %s/%s_%s.bcor", out_dir, aptamer, chr))
  cmd <- paste(cmd, sprintf("--bgen %s/%s_%s.bgen", out_dir, aptamer, chr))
  cmd <- paste(cmd, "--n-threads 1")
  system(command=cmd, wait=TRUE)

  # Convert LDstore bcor binary into usable table
  cmd <- "ldstore"
  cmd <- paste(cmd, sprintf("--bcor %s/%s_%s.bcor_1", out_dir, aptamer, chr))
  cmd <- paste(cmd, sprintf("--table %s/%s_%s.ld", out_dir, aptamer, chr))
  system(command=cmd, wait=TRUE)

  # Clean up
  system(sprintf("rm %s/%s_%s.bgen", out_dir, aptamer, chr), wait=TRUE)
  system(sprintf("rm %s/%s_%s.cis_window.bgen", out_dir, aptamer, chr), wait=TRUE)
  system(sprintf("rm %s/%s_%s.vars.txt", out_dir, aptamer, chr), wait=TRUE)
  system(sprintf("rm %s/%s_%s.bcor_1", out_dir, aptamer, chr), wait=TRUE)
}

# Combine per-chromosome scores into a single file
foreach(aptamer = unique(locmap$SOMAMER_ID)) %dopar% {
  chr_files <- list.files(path=out_dir, pattern=sprintf("%s_.*\\.ld", aptamer), full.names=TRUE)
  combined <- rbindlist(lapply(chr_files, fread))
  fwrite(combined, file=sprintf("%s/%s.ld", out_dir, aptamer), sep="\t", quote=FALSE)
  unlink(chr_files)
}

