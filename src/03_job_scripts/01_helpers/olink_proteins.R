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
olink_neu <- fread("data/INTERVAL/pre_qc_data/olink_proteomics/raw_data/high_dimensional_data/Olink_proteomics_neurology/gwasqc/olink_neu_qcgwas.csv", colClasses=c("aliquot_id"="character", "plate_id"="character"))

# Add IID to each panel
olink_inf[, join_key := aliquot_id]
olink_inf <- id_map[,.(IID=Affymetrix_gwasQC_bl, Olink_inf_gwasQC_24m)][olink_inf, on=.(Olink_inf_gwasQC_24m=join_key)]

olink_cvd2[, join_key := aliquot_id]
olink_cvd2 <- id_map[,.(IID=Affymetrix_gwasQC_bl, Olink_cvd2_gwasQC_24m)][olink_cvd2, on=.(Olink_cvd2_gwasQC_24m=join_key)]

olink_cvd3[, join_key := aliquot_id]
olink_cvd3 <- id_map[,.(IID=Affymetrix_gwasQC_bl, Olink_cvd3_gwasQC_24m)][olink_cvd3, on=.(Olink_cvd3_gwasQC_24m=join_key)]

olink_neu[, join_key := aliquot_id]
olink_neu <- id_map[,.(IID=Affymetrix_gwasQC_bl, Olink_neu_gwasQC_24m)][olink_neu, on=.(Olink_neu_gwasQC_24m=join_key)]

# Split off sample information
inf_sample_info <- olink_inf[, .(panel="inf", IID, aliquot_id, plate, Olink_inf_gwasQC_24m, Affymetrix_gwasQC_bl=IID)]
cvd2_sample_info <- olink_cvd2[, .(panel="cvd2", IID, aliquot_id, plate, Olink_cvd2_gwasQC_24m, Affymetrix_gwasQC_bl=IID)]
cvd3_sample_info <- olink_cvd3[, .(panel="cvd3", IID, aliquot_id, plate, Olink_cvd3_gwasQC_24m, Affymetrix_gwasQC_bl=IID)]
neu_sample_info <- olink_neu[, .(panel="neu", IID, aliquot_id, plate=plate_id, Olink_neu_gwasQC_24m, Affymetrix_gwasQC_bl=IID)]
sample_info <- rbind(inf_sample_info, cvd2_sample_info, cvd3_sample_info, neu_sample_info, fill=TRUE)

# Remove extraneous columns
olink_inf[, c("aliquot_id", "Olink_inf_gwasQC_24m", "plate") := NULL]
olink_cvd2[, c("aliquot_id", "Olink_cvd2_gwasQC_24m", "plate") := NULL]
olink_cvd3[, c("aliquot_id", "Olink_cvd3_gwasQC_24m", "plate") := NULL]
olink_neu[, c("aliquot_id", "Olink_neu_gwasQC_24m", "plate_id") := NULL]

# Filter to just people with IIDs
olink_inf <- olink_inf[!is.na(IID)]
olink_cvd2 <- olink_cvd2[!is.na(IID)]
olink_cvd3 <- olink_cvd3[!is.na(IID)]
olink_neu <- olink_neu[!is.na(IID)]

# Fix bad identifier
names(olink_inf)[71] <- "dner___q8nft8" # was "DNER___Q8NFT8Â", not sure what happend to this one

# Create combined table in long format
olink_inf_l <- melt(olink_inf, id.vars="IID", variable.name="protein", value.name="level")
olink_cvd2_l <- melt(olink_cvd2, id.vars="IID", variable.name="protein", value.name="level")
olink_cvd3_l <- melt(olink_cvd3, id.vars="IID", variable.name="protein", value.name="level")
olink_neu_l <- melt(olink_neu, id.vars="IID", variable.name="protein", value.name="level")
olink_all <- rbind(inf=olink_inf_l, cvd2=olink_cvd2_l, cvd3=olink_cvd3_l, neu=olink_neu_l, idcol="panel")

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

neu_info <- read.xlsx("data/INTERVAL/pre_qc_data/olink_proteomics/raw_data/high_dimensional_data/Olink_proteomics_neurology/rawData/20161008_Danesh_NPX_NEU.xlsx", sheet=1, colNames=FALSE)
neu_info <- as.data.table(t(neu_info[4:6,-1]))
neu_info <- neu_info[1:92]
setnames(neu_info, c("protein", "UniProt", "Olink_id"))
neu_info[, UniProt := gsub(",", ";", UniProt)] # make consistent with other panels
neu_info[, raw_id := paste0(gsub("OID00", "", Olink_id), "_", protein)] # make consistent with other panels
neu_info[, qc_id := tolower(paste0(gsub("\\.", "", make.names(protein)), "___", gsub("\\.", "", make.names(UniProt))))]

# For the neurology panel, we have to assign some qc_id's manually since they are non-standard format, or use
# gene symbols instead of protein symbols.
neu_info[UniProt == "Q08708", qc_id := "cd300c___q08708"] # not clm6___q08708
neu_info[UniProt == "P30533", qc_id := "lrpap1___p30533"] # not alpha2mrap___p30533
neu_info[UniProt == "Q92765", qc_id := "frzb___q92765"] # not sfrp3___q92765
neu_info[UniProt == "O00214", qc_id := "lgals8___o00214"] # not gal8___o00214
neu_info[UniProt == "P08473", qc_id := "mme___p08473"] # not nep___p08473
neu_info[UniProt == "O14793", qc_id := "mstn___o14793"] # not gdf8___o14793
neu_info[UniProt == "P56159", qc_id := "gfra1___p56159"] # not gfralpha1___p56159
neu_info[UniProt == "P15509", qc_id := "csf2ra___p15509"] # not gmcsfralpha___p15509
neu_info[UniProt == "P01138", qc_id := "ngf___p01138"] # not betangf___p01138
neu_info[UniProt == "P09919", qc_id := "csf3___p09919"] # not gcsf___p09919
neu_info[UniProt == "O60609", qc_id := "gfra3___o60609"] # not gdnfralpha3___o60609
neu_info[UniProt == "P37023", qc_id := "acvrl1___p37023"] # not skr3___p37023
neu_info[UniProt == "P78333", qc_id := "gpc5___p78333"] # not gcp5___p78333
neu_info[UniProt == "Q01344", qc_id := "il5ra___q01344"] # not il5ralpha___q01344
neu_info[UniProt == "P16234", qc_id := "pdgfra___p16234"] # not pdgfralpha___p16234
neu_info[UniProt == "P57087", qc_id := "jam2___p57087"] # not jamb___p57087
neu_info[UniProt == "Q9NR71", qc_id := "asah2___q9nr71"] # not ncdase___q9nr71
neu_info[UniProt == "Q9BZM5", qc_id := "ulbp2___q9bzm5"] # not n2dl2___q9bzm5
neu_info[UniProt == "Q8TDQ1", qc_id := "cd300lf___q8tdq1"] # not clm1___q8tdq1
neu_info[protein == "IL12", qc_id := "il12a_il12b___p29459_p29460"] # not il12___p29460p29459

# Make consistent column ordering
inf_info <- inf_info[, .(qc_id, protein, UniProt, raw_id, Olink_id)]
cvd2_info <- cvd2_info[, .(qc_id, protein, UniProt, raw_id, Olink_id)]
cvd3_info <- cvd3_info[, .(qc_id, protein, UniProt, raw_id, Olink_id)]
neu_info <- neu_info[, .(qc_id, protein, UniProt, raw_id, Olink_id)]

# make combined table
protein_info <- rbind(inf=inf_info, cvd2=cvd2_info, cvd3=cvd3_info, neu=neu_info, idcol="panel")

# Write out files
setnames(olink_all, c("protein", "level"), c("variable", "value"))
olink_all[, variable := paste0(panel, "_", variable)]
fwrite(olink_all[, .(IID, variable, value)], file=sprintf("%s/traits.tsv", out_dir), sep="\t", quote=FALSE)

protein_info[, variable := paste0(panel, "_", qc_id)]
fwrite(protein_info, file=sprintf("%s/trait_info.tsv", out_dir), sep="\t", quote=FALSE)


