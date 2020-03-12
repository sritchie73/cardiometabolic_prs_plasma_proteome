library(data.table)
library(openxlsx)

# make output directory
out_dir <- "analyses/processed_traits/olink_proteins"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Load data - come back and fix filepaths once data is migrated correctly
id_map <- fread("data/INTERVAL/project_1074/omicsMap.csv", colClasses="character", na.strings=c("NA", ""))
olink_inf <- fread("data/INTERVAL/pre_qc_data/olink_proteomics/raw_data/high_dimensional_data/Olink_proteomics_inf/gwasqc/olink_qcgwas_inf.csv", colClasses=c("aliquot_id"="character", "plate"="character"))
olink_cvd2 <- fread("data/INTERVAL/pre_qc_data/olink_proteomics/raw_data/high_dimensional_data/Olink_proteomics_cvd2/gwasqc/olink_qcgwas_cvd2.csv", colClasses=c("aliquot_id"="character", "plate"="character"))
olink_cvd3 <- fread("data/INTERVAL/pre_qc_data/olink_proteomics/raw_data/high_dimensional_data/Olink_proteomics_cvd3/gwasqc/olink_qcgwas_cvd3.csv", colClasses=c("aliquot_id"="character", "plate"="character"))

# Add IID to each panel
olink_inf[, join_key := aliquot_id]
olink_inf <- id_map[,.(IID=Affymetrix_gwasQC_bl, Olink_inf_gwasQC_24m)][olink_inf, on=.(Olink_inf_gwasQC_24m=join_key)]

olink_cvd2[, join_key := aliquot_id]
olink_cvd2 <- id_map[,.(IID=Affymetrix_gwasQC_bl, Olink_cvd2_gwasQC_24m)][olink_cvd2, on=.(Olink_cvd2_gwasQC_24m=join_key)]

olink_cvd3[, join_key := aliquot_id]
olink_cvd3 <- id_map[,.(IID=Affymetrix_gwasQC_bl, Olink_cvd3_gwasQC_24m)][olink_cvd3, on=.(Olink_cvd3_gwasQC_24m=join_key)]

# Split off sample information
inf_sample_info <- olink_inf[, .(panel="inf", IID, aliquot_id, plate, Olink_inf_gwasQC_24m, Affymetrix_gwasQC_bl=IID)]
cvd2_sample_info <- olink_cvd2[, .(panel="cvd2", IID, aliquot_id, plate, Olink_cvd2_gwasQC_24m, Affymetrix_gwasQC_bl=IID)]
cvd3_sample_info <- olink_cvd3[, .(panel="cvd3", IID, aliquot_id, plate, Olink_cvd3_gwasQC_24m, Affymetrix_gwasQC_bl=IID)]
sample_info <- rbind(inf_sample_info, cvd2_sample_info, cvd3_sample_info, neu_sample_info, fill=TRUE)

# Remove extraneous columns
olink_inf[, c("aliquot_id", "Olink_inf_gwasQC_24m", "plate") := NULL]
olink_cvd2[, c("aliquot_id", "Olink_cvd2_gwasQC_24m", "plate") := NULL]
olink_cvd3[, c("aliquot_id", "Olink_cvd3_gwasQC_24m", "plate") := NULL]

# Filter to just people with IIDs
olink_inf <- olink_inf[!is.na(IID)]
olink_cvd2 <- olink_cvd2[!is.na(IID)]
olink_cvd3 <- olink_cvd3[!is.na(IID)]

# Fix bad identifier
names(olink_inf)[71] <- "dner___q8nft8" # was "DNER___Q8NFT8Â", not sure what happend to this one

# Create combined table in long format
olink_inf_l <- melt(olink_inf, id.vars="IID", variable.name="protein", value.name="level")
olink_cvd2_l <- melt(olink_cvd2, id.vars="IID", variable.name="protein", value.name="level")
olink_cvd3_l <- melt(olink_cvd3, id.vars="IID", variable.name="protein", value.name="level")
olink_all <- rbind(inf=olink_inf_l, cvd2=olink_cvd2_l, cvd3=olink_cvd3_l, idcol="panel")

# Pull in raw data to get protein identifiers
inf_info <- read.xlsx("data/INTERVAL/pre_qc_data/olink_proteomics/raw_data/high_dimensional_data/Olink_proteomics_inf/rawData/20161007_Danesh_NPX_Inflammation.xlsx", sheet=1, colNames=FALSE)
inf_info <- as.data.table(t(inf_info[1:2,-1]))
inf_info <- inf_info[1:92]
setnames(inf_info, c("UniProt", "raw_id"))
inf_info[, UniProt := gsub("\\.", "", make.names(UniProt))] # remove irregular trailing whitespace
inf_info[, protein := gsub(".*_", "", raw_id)]
inf_info[, qc_id := tolower(paste0(gsub("\\.", "", make.names(protein)), "___", UniProt))]
inf_info[protein == "4E-BP1", qc_id := gsub("^x4", "", qc_id)] # fix since protein starts with number, which is missing in qced data identifier
inf_info[, Olink_id := paste0("OID00", gsub("_.*", "", raw_id))] # Add missing identifier present in neurology panel

cvd2_info <- read.xlsx("data/INTERVAL/pre_qc_data/olink_proteomics/raw_data/high_dimensional_data/Olink_proteomics_cvd2/rawData/20161007_Danesh_NPX_CVD2.xlsx", sheet=1, colNames=FALSE)
cvd2_info <- as.data.table(t(cvd2_info[1:2,-1]))
cvd2_info <- cvd2_info[1:92]
setnames(cvd2_info, c("UniProt", "raw_id"))
cvd2_info[, protein := gsub(".*_", "", raw_id)]
cvd2_info[, qc_id := tolower(paste0(gsub("\\.", "", make.names(protein)), "___", gsub("\\.", "", make.names(UniProt))))]
cvd2_info[, Olink_id := paste0("OID00", gsub("_.*", "", raw_id))]  # Add missing identifier present in neurology panel

cvd3_info <- read.xlsx("data/INTERVAL/pre_qc_data/olink_proteomics/raw_data/high_dimensional_data/Olink_proteomics_cvd3/rawData/20161007_Danesh_NPX_CVD3_update.xlsx", sheet=1, colNames=FALSE)
cvd3_info <- as.data.table(t(cvd3_info[1:2,-1]))
cvd3_info <- cvd3_info[1:92]
setnames(cvd3_info, c("UniProt", "raw_id"))
cvd3_info[, protein := gsub(".*_", "", raw_id)]
cvd3_info[, qc_id := tolower(paste0(gsub("\\.", "", make.names(protein)), "___", gsub("\\.", "", make.names(UniProt))))]
cvd3_info[, Olink_id := paste0("OID00", gsub("_.*", "", raw_id))]  # Add missing identifier present in neurology panel

# Make consistent column ordering
inf_info <- inf_info[, .(qc_id, protein, UniProt, raw_id, Olink_id)]
cvd2_info <- cvd2_info[, .(qc_id, protein, UniProt, raw_id, Olink_id)]
cvd3_info <- cvd3_info[, .(qc_id, protein, UniProt, raw_id, Olink_id)]

# make combined table
protein_info <- rbind(inf=inf_info, cvd2=cvd2_info, cvd3=cvd3_info, idcol="panel")

# Write out files
setnames(olink_all, c("protein", "level"), c("variable", "value"))
olink_all[, variable := paste0(panel, "_", variable)]
fwrite(olink_all[, .(IID, variable, value)], file=sprintf("%s/traits.tsv", out_dir), sep="\t", quote=FALSE)

protein_info[, variable := paste0(panel, "_", qc_id)]
fwrite(protein_info, file=sprintf("%s/trait_info.tsv", out_dir), sep="\t", quote=FALSE)


