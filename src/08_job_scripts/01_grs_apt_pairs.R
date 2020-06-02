# Get a list of GRS-aptamer pairs that:
#  (1) Are significantly association (FDR < 0.05)
#  (2) Where the aptamer's pQTLs attenuate the GRS-association
library(data.table)
library(foreach)
source("src/utilities/prot_pval.R")

out_dir <- "analyses/grs_pqtl_removed"

# Get GRSs to test from stdin
args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() != 0L) {  # if __name__ == '__main__':
  options(datatable.print.nrows=30)
  options(width=170)
  stopifnot(exists("view_file"))
  grs_info <- fread(view_file)
  GRSs <- grs_info[, GRS_name]
  view <- gsub(".txt$", "", basename(view_file))
  cat("Comparing", view, "GRSs.\n")
} else if (length(args) == 0) {
  GRSs <- list.files("analyses/GRS_profiles/")
  view <- "all"
  cat("Comparing all GRSs.\n")
} else if (length(args) == 1) {
  view_file <- args[1]
  grs_info <- fread(view_file)
  GRSs <- grs_info[, GRS_name]
  view <- gsub(".txt$", "", basename(view_file))
  cat("Comparing", view, "GRSs.\n")
} else {
  stop("Unexpected extra arguments.")
}

# Load associations
assocs <- foreach(grs = GRSs) %do% {
  fread(sprintf("analyses/univariate_associations/%s/somalogic_proteins/associations.tsv", grs))
}
names(assocs) <- GRSs
assocs <- rbindlist(assocs, idcol="grs")

# Load protein information to filter to "good" aptamers
prot_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
prot_info <- prot_info[Type == "Protein"]

# Add information, filter, average associations across aptamers, and filter at FDR < 0.05:
assocs <- assocs[trait %in% prot_info$variable]
assocs[prot_info, on = .(trait=variable), Target := Target]
assocs <- assocs[, .(pval = prot_pvalue(pval, beta)), by=.(grs, Target)]
assocs[, fdr := p.adjust(pval, method="fdr"), by=grs]
assocs <- assocs[fdr < 0.1, .(grs, Target)]

# Add back aptamer set
assocs <- assocs[prot_info[,.(Target, variable)], on = .(Target), nomatch=0]

# Output pairs list
fwrite(assocs[,.(grs, variable)], quote=FALSE, sep="\t", col.names=FALSE, 
       file=sprintf("%s/pairs_to_run.txt", out_dir))
