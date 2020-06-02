# Output gene lists for gene set enrichment analysis with gProfiler
library(data.table)
library(foreach)

out_dir <- "analyses/pub/cardiometabolic_proteins/review1"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Load in protein information
info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
# Filter to "good" aptamers
info <- info[Type == "Protein"]
# Get list of unique genes to use as the background set
background <- info[, .(Gene=unlist(strsplit(Gene.Name, "\\|")))]
background <- background[, .(Gene = unlist(strsplit(Gene, ",")))]
background <- unique(background)
fwrite(background, quote=FALSE, col.names=FALSE, file=sprintf("%s/background_set.txt", out_dir))

# Load in FDR significant associations to run enrichment analysis for
view_file <- "views/cardiometabolic.txt"
grs_info <- fread(view_file)
GRSs <- grs_info[, GRS_name]

assocs <- foreach(grs = GRSs) %do% {
  fread(sprintf("analyses/univariate_associations/%s/somalogic_proteins/associations.tsv", grs))
}
names(assocs) <- GRSs
assocs <- rbindlist(assocs, idcol="grs")

# Add information, filter, average associations across aptamers, and filter at FDR < 0.05:
assocs <- assocs[trait %in% info$variable]
assocs <- assocs[info, on = .(trait=variable), nomatch=0]
prot_assocs <- assocs[, .(Beta = mean(beta), L95 = mean(l95), U95 = mean(u95), Pvalue = mean(pval)),
                        by = .(PRS=grs, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, Gene=Gene.Name)]
prot_assocs[, FDR := p.adjust(Pvalue, method="fdr"), by=PRS]
prot_assocs <- prot_assocs[FDR < 0.05]
prot_assocs[Gene == "PDE4D", c("Gene", "UniProt") := .("PDE4D|PDE4A", "Q08499|P27815")]

# Split out entries where there are multiple proteins being measured by the same ap
for (grs in GRSs) {
  this_assocs <- prot_assocs[PRS == grs]
  if (nrow(this_assocs) > 0) {
    sig_genes <- this_assocs[,.(Gene = unlist(strsplit(Gene, "\\|")))]
    sig_genes <- sig_genes[, .(Gene = unlist(strsplit(Gene, ",")))]
    sig_genes <- unique(sig_genes)
    fwrite(sig_genes, quote=FALSE, col.names=FALSE, file=sprintf("%s/%s_sig_set.txt", out_dir, grs))
  } 
}

