library(data.table)
library(openxlsx)
source("src/utilities/format_pval.R")

# Load all sets of associations
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
hes_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")
mr_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/collated_mr.txt")

# Extract protein level associations
prs_assocs <- unique(prs_assocs[,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])
hes_assocs <- unique(hes_assocs[,.(endpoint, Target, UniProt, Gene, HR, HR.L95, HR.U95, HR.P, HR.FDR)])
mr_assocs <- unique(mr_assocs[,.(GWAS, Target, UniProt, Gene, MR=exp(mr_estimate), MR.L95=exp(mr_L95), MR.U95=exp(mr_U95), MR.P=mr_pval, 
                                 MR.FDR=mr_fdr, MR.pleiotropy=pleiotropy_pval, colocalization)])

# Map between outcomes
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  endpoint = c("Atrial fibrillation", "Myocardial infarction", "End stage renal disease", "Ischaemic stroke", "Diabetes"),
  GWAS = c("Afib", "CAD", "CKD", "StrokeIS", "T2DadjBMI")
)

prs_assocs[map, on = .(PRS), c("GWAS", "endpoint") := .(GWAS, endpoint)]
hes_assocs[map, on = .(endpoint), c("PRS", "GWAS") := .(PRS, GWAS)]
mr_assocs[map, on = .(GWAS), c("PRS", "endpoint") := .(PRS, endpoint)]

# Filter to FDR < 0.05 PRS to protein associations
prs_assocs <- prs_assocs[PRS.FDR < 0.05]

# Build comparison table
comp <- merge(prs_assocs, hes_assocs, by=c("PRS", "GWAS", "endpoint", "Target", "UniProt", "Gene"), all.x=TRUE)
comp <- merge(comp, mr_assocs, by=c("PRS", "GWAS", "endpoint", "Target", "UniProt", "Gene"), all.x=TRUE)

format_beta <- function(n) {
  sprintf("%s", ifelse(abs(n) > 0.1, round(n, digits=2), round(n, digits=3)))
}

format_hr <- function(n) {
  sprintf("%s", ifelse(abs(n-1) > 0.1, round(n, digits=2), round(n, digits=3)))
}

# Convert to text
comp[, PRS.Beta := format_beta(PRS.Beta)]
comp[, PRS.95CI := sprintf("[%s, %s]", format_beta(PRS.L95), format_beta(PRS.U95))]
comp[, PRS.P := my_format_pval(PRS.P)]
comp[, PRS.FDR := my_format_pval(PRS.FDR)]

comp[, HR := format_hr(HR)]
comp[, HR.95CI := sprintf("[%s, %s]", format_hr(HR.L95), format_hr(HR.U95))]
comp[, HR.P := my_format_pval(HR.P)]
comp[, HR.FDR := my_format_pval(HR.FDR)]

comp[, MR := format_hr(MR)]
comp[, MR.95CI := sprintf("[%s, %s]", format_hr(MR.L95), format_hr(MR.U95))]
comp[, MR.P := my_format_pval(MR.P)]
comp[, MR.FDR := my_format_pval(MR.FDR)]
comp[, MR.pleiotropy := my_format_pval(MR.pleiotropy)]
comp[, colocalization := ifelse(colocalization, "Yes", "No")]

comp[HR == "NA", c("HR", "HR.95CI", "HR.P", "HR.FDR") := ""]
comp[MR == "NA", c("MR", "MR.95CI", "MR.P", "MR.FDR", "MR.pleiotropy", "colocalization") := ""]

# Flag pQTL status
prot_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
prot_info <- prot_info[Type == "Protein"]
prot_info <- prot_info[, .(SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, Gene=Gene.Name)]

sun <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=5, startRow=3)
sun <- as.data.table(sun)
head_row <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, rows=5:6, fillMergedCells=TRUE)
head_row <- gsub(".NA$", "", paste(colnames(head_row), as.vector(head_row[1,]), sep="."))
pQTL_info <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, startRow=6)
colnames(pQTL_info) <- head_row
pQTL_info <- as.data.table(pQTL_info)
sun[pQTL_info, on = .(SOMAmer.ID, Sentinel.variant=`Sentinel.variant*`), type := `cis/.trans`]
sun <- sun[,.(SOMAMER_ID=SOMAmer.ID, rsid=Conditional.variant, type)]
sun <- sun[SOMAMER_ID %in% prot_info$SOMAMER_ID]
sun <- sun[prot_info, on = .(SOMAMER_ID), nomatch=0]

new_trans <- fread("analyses/pub/cardiometabolic_proteins/review2/new_trans_pQTLs.txt")
new_trans <- new_trans[Aptamer %in% prot_info$Aptamer]
new_trans <- new_trans[,.SD[which.min(P)],by=.(Target, UniProt, Gene, Aptamer)]

cis <- fread("analyses/pub/cardiometabolic_proteins/review2/cis_pQTLs.txt")
cis <- cis[Aptamer %in% prot_info$Aptamer]

prot_info[, sun_trans := ifelse(SOMAMER_ID %in% sun[type == "trans", SOMAMER_ID], TRUE, FALSE)]
prot_info[, sun_cis := ifelse(SOMAMER_ID %in% sun[type == "cis", SOMAMER_ID], TRUE, FALSE)]
prot_info[, new_trans := ifelse(Aptamer %in% new_trans$Aptamer, TRUE, FALSE)]
prot_info[, new_cis := ifelse(Aptamer %in% cis$Aptamer, TRUE, FALSE)]

comp[, type := "(0) No pQTLs"]
comp[Target %in% prot_info[(sun_trans) | (new_trans), Target], type := "(1) Trans pQTLs only"]
comp[Target %in% prot_info[(sun_cis) | (new_cis), Target], type := "(2) At least one cis-pQTL"]
comp[MR != "", type := "(3) Three or more cis-pQTLs"]

# Add polygenicity information
polygenicity <- fread("analyses/pub/cardiometabolic_proteins/review2/polygenicity.txt")
comp[polygenicity, on = .(PRS, Target), polygenicity := paste0(round(i.pct_removed*100), "%")]
comp[polygenicity, on = .(PRS, Target), LD_blocks := LD_blocks_removed]

# Arrange columns
comp <- comp[, .(Target, UniProt, Gene, PRS, PRS.Beta, PRS.95CI, PRS.P, PRS.FDR, polygenicity,
                 endpoint, HR, HR.95CI, HR.P, HR.FDR,
                 type, GWAS, MR, MR.95CI, MR.P, MR.FDR, MR.pleiotropy, colocalization)]
comp[polygenicity == "0%", c("PolyMR", "PolyMR.95CI", "PolyMR.P", "PolyMR.pleiotropy") := ""]

fwrite(comp, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/prs_prot_compare.txt")
