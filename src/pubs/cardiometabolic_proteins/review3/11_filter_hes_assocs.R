library(data.table)

prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[FDR < 0.05,.(PRS, Gene)])
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  endpoint = c("Atrial fibrillation", "Myocardial infarction", "End stage renal disease", "Ischaemic stroke", "Diabetes")
)
prs_assocs[map, on = .(PRS), endpoint := i.endpoint]


hes_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")
hes_assocs[map, on = .(endpoint), PRS := i.PRS]
hes_assocs <- hes_assocs[prs_assocs, on = .(PRS, endpoint, Gene)]
hes_assocs <- hes_assocs[!(PRS %in% c("CKD_PRS", "IS_PRS"))]

fwrite(hes_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/filtered_hes_assocs.txt")

