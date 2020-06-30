library(data.table)
library(RNOmni)

# make output directory
out_dir <- "analyses/processed_traits/olink_ngs"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# load in new olink NGS data:
ngs <- fread("data/INTERVAL/olink_ngs/20200292_Danesh_NPX_2020-06-03.csv")

# Get plate and well:
ngs[, well := gsub(".*_", "", SampleID)]
ngs[, plate := gsub("_.*", "", SampleID)]

# Drop samples in column 1 of each plate - apparently these samples had problems
ngs <- ngs[!grepl("01", well)]

# Map samples IDs to affymetrix IDs and filter to samples with genotype data:
idmap1 <- fread("data/INTERVAL/project_1074/2020-06-19/Olink_NGS.csv")
idmap2 <- fread("data/INTERVAL/project_1074/2020-06-19/INTERVAL_OmicsMap_20200619.csv")
idmap2[idmap1, on = .(identifier), Olink_ngs_24m := OLINK_NGS]

# Give the new olink data the affymetrix IDs, also removing the 2 samples that overlap with the somalogic data:
ngs[idmap2[!is.na(Affymetrix_gwasQC_bl) & is.na(soma4000_gwasQC_bl)], on = .(SampleID=Olink_ngs_24m), IID := Affymetrix_gwasQC_bl]
ngs <- ngs[!is.na(IID)]

# Split out into information and protein measurements:
prot_info <- ngs[,.(Panel, Panel_Version, OlinkID, UniProt, Gene=Assay, MissingFreq, variable=paste0(Panel, "_", UniProt))]
prot_info <- unique(prot_info)

sample_info <- ngs[!is.na(IID), .(IID, plate, well, SampleID)]
sample_info <- unique(sample_info)

raw <- ngs[,.(IID, variable = paste0(Panel, "_", UniProt), value=NPX)]
tags <- ngs[,.(IID, variable = paste0(Panel, "_", UniProt), QC_Warning, LOD)]

# Load phenotype information and derive covariates to adjust for
pheno <- fread("data/INTERVAL/project_1074/INTERVALdata_17JUL2018.csv")
idmap3 <- fread("data/INTERVAL/project_1074/omicsMap.csv")
pheno[idmap3, on = .(identifier), IID := Affymetrix_gwasQC_bl]

pheno[, daysToProcess := as.IDate(processDate_24m, format="%d%b%Y") - as.IDate(attendanceDate_24m, format="%d%b%Y")]
pheno[, followup_24m := as.IDate(attendanceDate_24m, format="%d%b%Y") - as.IDate(attendanceDate, format="%d%b%Y")]
pheno[, agePulse_24m := agePulse + round(followup_24m/365*10)/10]

# Adjust for age, sex, days between blood draw and sample processing, and plate
raw[pheno, on = .(IID), c("age", "sex", "daysToProcess") := .(agePulse_24m, sexPulse, daysToProcess)]
raw[sample_info, on = .(IID), plate := plate]
raw <- raw[!is.na(value)]
raw[, adj := lm(value ~ scale(age) + factor(sex) + daysToProcess + factor(plate))$residuals, by=.(variable)]

# Inverse rank normalise:
raw[, adj := rankNorm(adj), by = .(variable)]

# Write out:
fwrite(raw[,.(IID, variable, value=adj)], sep="\t", quote=FALSE, file=sprintf("%s/traits.tsv", out_dir))
fwrite(tags, sep="\t", quote=FALSE, file=sprintf("%s/qc_tags.tsv", out_dir))
fwrite(prot_info, sep="\t", quote=FALSE, file=sprintf("%s/trait_info.tsv", out_dir))
fwrite(sample_info, sep="\t", quote=FALSE, file=sprintf("%s/sample_info.tsv", out_dir))

