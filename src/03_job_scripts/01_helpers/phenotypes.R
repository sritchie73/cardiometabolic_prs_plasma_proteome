library(data.table)

# Not in its own directory so it doesnt get picked up by the polygenic association scan scripts
out_dir <- "analyses/processed_traits"

# Load phenotype data, set identifier to the genetic identifier, and remove blood cell traits
pheno <- fread("data/INTERVAL/project_1074/INTERVALdata_17JUL2018.csv", colClasses=c("identifier"="IID"))
id_map <- fread("data/INTERVAL/project_1074/omicsMap.csv", colClasses="character", na.strings=c("NA", ""))

pheno[id_map, on = .(identifier), IID := Affymetrix_gwasQC_bl]
pheno <- pheno[, c("IID", names(pheno)[1:57], names(pheno)[147:149]), with=FALSE]

fwrite(pheno, file=sprintf("%s/phenotypes.tsv", out_dir), sep="\t", quote=FALSE)

