library(data.table)

sens <- rbind(idcol="type", fill=TRUE, use.names=TRUE,
  "Adjusting for circardian effects" = fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_time_assocs.txt"),
  "Adjusting for seasonal effects" = fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_date_assocs.txt"),
  "Including prevalent cardiometabolic disease" = fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs_with_prev.txt"),
  "Adjusting for BMI" = fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_bmi_assocs.txt")
)
sens[, type := factor(type, levels=unique(type))]
sens <- sens[is.na(Coef) | Coef == "PRS"]

prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[FDR < 0.05,.(PRS, Target, UniProt, Gene)])
sens <- sens[prs_assocs, on = .(PRS, Target, UniProt, Gene)]

fwrite(sens, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/filtered_sensitivity.txt")


