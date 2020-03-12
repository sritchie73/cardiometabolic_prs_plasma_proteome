library(data.table)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
grs_name <- args[1]
aptamer <- args[2]
grs_resource_dir <- file.path("data/GRS_resources", grs_name)
analysis_dir <- file.path("analyses/grs_pqtl_removed", grs_name, aptamer)

# ----------------------------------------------------------------------------
# First build a table containing the number of variants passing QC in each GRS
# ----------------------------------------------------------------------------

# Get the number of lines in the files for each GRS.
GRS_variants <- as.numeric(system(paste("zcat", file.path(grs_resource_dir, "grs_weights.txt.gz"), "| wc -l"), intern=T))
GRS_variants <- GRS_variants - 1

# And the number of variants passing QC for each chromosome
GRS_passing_qc <- system(sprintf("wc -l %s/*.sscore.vars | tail -1", analysis_dir), intern=TRUE)
GRS_passing_qc <- gsub(" total", "", GRS_passing_qc)
GRS_passing_qc <- as.numeric(GRS_passing_qc)

msg <- paste0(format(GRS_passing_qc, big.mark=","), " (", round(GRS_passing_qc/GRS_variants*100), "%)",
              " of the ", format(GRS_variants, big.mark=","),
              " variants from the GRS passed QC when calculating each person's score.")
cat(msg, "\n", file=file.path(analysis_dir, "score_summary.txt"))

# -------------------------------------------------------------------
# For each GRS, combine the per-chromosome scores for each individual
# -------------------------------------------------------------------

# Load the per-chromosome score files
chr_files <- list.files(".sscore(.gz)?$", path=analysis_dir, full.names=TRUE)
chr_scores <- lapply(chr_files, fread)
chr_scores <- rbindlist(chr_scores)
setnames(chr_scores, "#IID", "IID")

# Each sample has a per-chromosome score, SCORE1_SUM, comprising the sum of the 
# meta-GRS coefficients * effect allele presence in that chromosome. To get the 
# final per individual score, we sum those scores across all chromosome, then
# standardise across the population
chr_scores <- chr_scores[, .(non_missing_alleles=sum(NMISS_ALLELE_CT), score_sum=sum(SCORE1_SUM)), by = IID]
chr_scores[, score := as.vector(scale(score_sum))] # standardise the GRS.
chr_scores[, missing_alleles := (GRS_variants * 2) - non_missing_alleles]

# write out per-individual scores to file:
fwrite(chr_scores, file.path(analysis_dir, "profile.sscore"), sep="\t")
