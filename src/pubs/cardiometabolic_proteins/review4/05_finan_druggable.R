library(data.table)
library(openxlsx)

# Get list of PGS-associated proteins with evidence for causal mediation
pgs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review3/compare_estimates.txt")
pgs_assocs <- unique(pgs_assocs[Causal.P < 0.05,.(Gene)])

# Load druggable targets
druggable <- read.xlsx("data/Finan_et_al_2017/Table S1.xlsx", sheet="Data")
setDT(druggable)

# Intersect
pgs_druggable <- druggable[pgs_assocs, on = .(hgnc_names=Gene), nomatch=0]
pgs_druggable <- pgs_druggable[order(hgnc_names)][order(druggability_tier)]

fwrite(pgs_druggable, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/druggable_targets.txt")
