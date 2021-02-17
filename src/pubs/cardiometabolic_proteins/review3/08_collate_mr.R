library(data.table)

# Load data
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
mr_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/collated_mr.txt")

# Extract protein level associations
prs_assocs <- unique(prs_assocs[,.(PRS, Gene, PRS.FDR=FDR)])
mr_assocs <- unique(mr_assocs[,.(GWAS, Gene, MR=exp(mr_estimate), MR.L95=exp(mr_L95), MR.U95=exp(mr_U95), MR.P=mr_pval,
                                 MR.pleiotropy=pleiotropy_pval, colocalization)])

# Map between outcomes
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  GWAS = c("Afib", "CAD", "CKD", "StrokeIS", "T2DadjBMI")
)

prs_assocs[map, on = .(PRS), GWAS := i.GWAS]
mr_assocs[map, on = .(GWAS), PRS := i.PRS]

# Filter to FDR < 0.05 PRS to protein associations
prs_assocs <- prs_assocs[PRS.FDR < 0.05, .(PRS, GWAS, Gene)]
mr_assocs <- mr_assocs[prs_assocs, on = .(PRS, GWAS, Gene), nomatch=0]

# Get estimates from specific MR methods
mr2 <- fread("analyses/pub/cardiometabolic_proteins/review2/collated_mr.txt")
mr2 <- unique(apt_mr[,.(GWAS, Gene, Method, MR=exp(mr_estimate.protein), MR.L95=exp(mr_L95.protein), MR.U95=exp(mr_U95.protein), MR.P=mr_pval.protein)])
mr_assocs <- mr_assocs[mr2, on = .(GWAS, Gene), nomatch=0]

# Get estimates for specific aptamers
mr3 <- fread("analyses/pub/cardiometabolic_proteins/review2/collated_mr.txt")
mr3 <- apt_mr[,.(GWAS, Gene, Method, Aptamer, MR=exp(mr_estimate.aptamer), MR.L95=exp(mr_L95.aptamer), MR.U95=exp(mr_U95.aptamer), MR.P=mr_pval.aptamer)]
mr_assocs <- mr_assocs[mr3, on = .(GWAS, Gene, Method), nomatch=0]

# Get row order
mro <- unique(mr_assocs[,.(PRS, GWAS, Gene, MR.P)])
mro <- mro[order(MR.P)]
mr_assocs <- mr_assocs[mro[,.(PRS, GWAS, Gene)], on = .(PRS, GWAS, Gene)]

# write out
fwrite(mr_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/collated_mr.txt")

# Load and collate instruments
ivs <- fread("analyses/pub/cardiometabolic_proteins/review2/collated_ivs.txt")
ivs[map, on = .(GWAS), PRS := i.PRS]
ivs <- ivs[mro[,.(PRS, GWAS, Gene)], on = .(PRS, GWAS, Gene), nomatch=0]
fwrite(ivs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/collated_ivs.txt")


