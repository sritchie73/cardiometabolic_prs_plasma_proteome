library(data.table)
library(openxlsx)

pQTLs <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=5, startRow=3)
pQTLs <- as.data.table(pQTLs)
pQTLs[, variable := tolower(gsub("\\.", "", SOMAmer.ID))]
pQTLs <- pQTLs[,.(variable, variant=Conditional.variant)]
fwrite(pQTLs, file="analyses/processed_traits/somalogic_proteins/qtl_list.txt", sep="\t", quote=FALSE)
fwrite(pQTLs[,.(unique(variant))],
       file="analyses/processed_traits/somalogic_proteins/qtl_variants.txt", 
       quote=FALSE, col.names=FALSE)
