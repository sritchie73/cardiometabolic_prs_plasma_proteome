library(data.table)
library(foreach)
library(openxlsx)

# Load pair file
pairs <- fread("analyses/grs_pqtl_removed/pairs_to_run.txt", header=FALSE)
setnames(pairs, c("PRS", "variable"))

# Load all files
info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
pQTLs <- as.data.table(read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=5, startRow=3))
bim <- foreach(chrIdx = 1:22, .combine=rbind) %do% {
  bim <- fread(sprintf("analyses/grs_pqtl_removed/chr%s.bim", chrIdx), header=FALSE)
  bim[, .(chr=V1, var_id=V2, pos=V4)]
}

# Filter protein info and pQTLs to the aptamers of interest
info <- info[variable %in% pairs$variable]
pQTLs[info, on = .(SOMAmer.ID = SOMAMER_ID), variable := variable]
pQTLs <- pQTLs[!is.na(variable)]

# For the protein info sheet, split out genomic locations if multiple
# e.g. due to multiple genes
info <- info[, .(
  chr=strsplit(chr, ",")[[1]],
  start=as.numeric(strsplit(start, ",")[[1]])
  ), by=variable]

# Filter to autosomes
info[, chr := suppressWarnings(as.numeric(chr))]
info <- info[!is.na(chr)]

# Determine 1 MB window around each gene:
info[, window_start := start - 1e6]
info[, window_end := start + 1e6]

# Determine 1 MB window around each pQTL
pQTLs[, window_start := Conditional.variant.Pos - 1e6]
pQTLs[, window_end := Conditional.variant.Pos + 1e6]

# combine exclusion regions:
excl_reg <- rbind(
  info[, .(variable, chr, window_start, window_end)],
  pQTLs[, .(variable, chr=Conditional.variant.Chr, window_start, window_end)])

# Define high complexity regions
# See flashpca exclusion regions: https://github.com/gabraham/flashpca
# coordinates are HG29
complex_ld <- data.table(
  region_chr=c(5, 6, 8, 11),
  region_start=c(44000000, 25000000, 8000000, 45000000),
  region_end=c(51500000, 33500000, 12000000, 57000000),
  region_name=c("r1", "MHC", "r3", "r4")
)

# If 1 MB window inside any region, redefine window as region
excl_reg[complex_ld, on = .(chr = region_chr, window_start > region_start, window_end < region_end),
          c("window_start", "window_end") := .(region_start, region_end)]
# If 1 MB window overlaps the region, extend the window:
excl_reg[complex_ld, on = .(chr = region_chr, window_start < region_end, window_end > region_end),
          c("window_start") := .(region_start)]
excl_reg[complex_ld, on = .(chr = region_chr, window_end > region_start, window_start < region_start),
          c("window_end") := .(region_end)]

# Extract variants to exclude for each test from the bim table
tags <- bim[excl_reg, on = .(chr, pos > window_start, pos < window_end), nomatch=0,
            .(variable, chr, var_id, pos=x.pos)]
tags <- unique(tags) # duplicates may arise where windows overlap


# Write out list of cis-variants to exclude for each pair of tests
foreach(rIdx = pairs[,.I]) %do% {
  PRS <- pairs[rIdx, PRS]
  var <- pairs[rIdx, variable]
  out_dir <- sprintf("analyses/grs_pqtl_removed/%s/%s/", PRS, var)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  fwrite(tags[variable == var, .(var_id)], col.names=FALSE, quote=FALSE, file=sprintf("%s/pqtl_tags.txt", out_dir))
  fwrite(tags[variable == var, .(unique(chr))], col.names=FALSE, quote=FALSE, file=sprintf("%s/chrs_to_filter.txt", out_dir))
}
