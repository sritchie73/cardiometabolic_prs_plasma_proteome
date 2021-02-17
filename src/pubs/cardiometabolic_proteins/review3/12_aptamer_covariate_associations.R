library(data.table)
library(foreach)
library(doMC)

# This script summarises the distribution of each aptamer and summarises their 
# associations with each covariate.
out_dir <- "analyses/pub/cardiometabolic_proteins/review3"

# set up parallel environment
registerDoMC(as.integer(Sys.getenv("SLURM_CPUS_ON_NODE")))

# Load in raw somalogic data (2 batches)
batch1 <- fread("data/INTERVAL/pre_qc_data/somalogic_proteomics/raw_data/SOMAS_r1_Raw_data_for_adam.tsv")
batch2 <- fread("data/INTERVAL/pre_qc_data/somalogic_proteomics/raw_data/SOMAS_r2_Raw_data_for_adam.tsv")

# To long format
batch1 <- melt(batch1, id.vars=names(batch1)[1:25], variable.name="Aptamer.ID")
batch2 <- melt(batch2, id.vars=names(batch2)[1:26], variable.name="Aptamer.ID")
raw <- rbind(batch1, batch2, idcol="batch", fill=TRUE)

# Load in post-QC data used for the rest of the study (to get samples passing QC)
gwasqc <- fread("data/INTERVAL/post_qc_data/phenotype/somalogic_proteomics/gwasqc/somalogic_qcgwas_4000.csv")
gwasqc <- gwasqc[, .(aliquot_id, samplegroup, batch)]

# Load in ID mapping files
omics_idmap <- fread("data/INTERVAL/project_1074/omicsMap.csv")
somas_idmap <- fread("data/SunB_etal_2018/sampleids/SOMAS_Samples_allids_plinkmapped.csv")

# Filter raw data to samples present in the post-QC dataset and also in the project id mapping file
somas_idmap <- unique(somas_idmap[,.(SOMAscan_ALIQUOT, PlasmaID)]) # some samples have longitudinal repeats
omics_idmap <- omics_idmap[!is.na(Affymetrix_gwasQC_bl) & !is.na(soma4000_gwasQC_bl)]
somas_idmap <- somas_idmap[PlasmaID %in% omics_idmap$soma4000_gwasQC_bl]
gwasqc <- gwasqc[somas_idmap, on = .(aliquot_id=PlasmaID)]
raw <- raw[gwasqc, on = .(SampleId=SOMAscan_ALIQUOT, batch)]

# Some aptamers have different IDs in the two batches
soma_info <- fread("data/INTERVAL/post_qc_data/phenotype/somalogic_proteomics/Documents/Somamer_info.tsv")
raw[, Aptamer.ID := as.vector(Aptamer.ID)] # otherwise factor
raw[, Seq.ID := sapply(strsplit(Aptamer.ID, split="\\."), function(x) { paste(x[(length(x)-2):(length(x)-1)], collapse="-") })]
raw[soma_info, on = .(Seq.ID=Seq), SOMAMER_ID := i.SOMAMER_ID]

# Extract out just a data.table of the raw protein levels, using the genotype data identifier so we 
# can map to the other information:
raw <- raw[omics_idmap, on = .(aliquot_id=soma4000_gwasQC_bl), .(IID=Affymetrix_gwasQC_bl, Seq.ID, value)]

# Filter to just the aptamers analysed:
prot_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
prot_info <- prot_info[Type == "Protein"]
raw <- raw[Seq.ID %in% prot_info$SeqId]

# Filter and fix the protein information
prot_info <- prot_info[,.(SeqId, Gene=Gene.Name, Protein=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot)]

prot_info[Protein == "14-3-3 protein family", UniProt := "P61981|Q04917"]
prot_info[Protein == "Induced myeloid leukemia cell differentiation protein Mcl-1", c("UniProt", "Gene") := .("Q07820", "MCL1")]
prot_info[Protein == "Protein delta homolog 1", c("UniProt", "Gene") :=  .("P80370", "DLK1")]
prot_info[Protein == "Stromal cell-derived factor 1", c("UniProt", "Gene") := .("P48061", "CXCL12")]

# Filter to samples without prevalent cardiometabolic disease
source("src/pubs/cardiometabolic_proteins/review2/load_HES.R")
hes <- hes[code_type == "primary+secondary"]
source("src/pubs/cardiometabolic_proteins/review2/cardiometabolic_events.R")
prev <- unique(hes[prevalent == 1 & phenotype %in% cardiometabolic, .(IID)])
raw <- raw[!prev, on = .(IID)]
raw <- raw[IID %in% unique(hes$IID)] # drops 1 sample with no EHR linkage

# Load phenotype data and build derived variables
pheno <- fread("analyses/processed_traits/phenotypes.tsv", na.strings=c("NA", ""))
pheno <- pheno[IID %in% unique(raw$IID)]

pheno[wt_bl == 777, wt_bl := NA] # bad coding
pheno[, bmi_raw := wt_bl/ht_bl^2]
pheno[, bmi := bmi_raw]
pheno[ht_bl < 1.47, bmi := NA_real_] # clinical cutoff for dwarfism
pheno[ht_bl > 2.1, bmi := NA_real_] # clinical cutoff for gigantism
pheno[wt_bl < 50 | wt_bl > 160, bmi := NA_real_] # NHS restrictions for weight

# Get the season (10-date buckets) the sample was taken
pheno[, date := as.POSIXct(attendanceDate, format="%d%b%Y")]
pheno[, year_day := yday(date)]
pheno[, year := year(date)]
pheno[year == 2012, from_2012 := year_day]
pheno[year == 2013, from_2012 := 366L + year_day]
pheno[year == 2014, from_2012 := 366L + 365L + year_day]
pheno[, from_start := from_2012 - min(from_2012) + 1L]
pheno[, date_bin := ceiling(from_start/(max(from_start)/10))]

# Get information about date bins:
date_bins <- pheno[,.N,by=date_bin]
date_start <- pheno[order(date), .SD[1], by=date_bin][,.(date_bin, start=date, start_n=from_2012)]
date_bins <- date_bins[date_start, on=.(date_bin), nomatch=0]
date_end <- date_start[, .(date_bin=date_bin-1, end=start, end_n=start_n)]
date_last <- pheno[order(date)][.N][, .(date_bin, end=date, end_n=from_2012)]
date_end <- date_end <- rbind(date_end[-1], date_last)
date_bins <- date_bins[date_end, on=.(date_bin), nomatch=0]
date_bins[, duration := end_n - start_n]
date_bins[, c("start_n", "end_n") := NULL]
date_bins <- date_bins[order(date_bin)]

# Code factor levels
reference_bin <- date_bins[which.max(N), date_bin]
pheno[, date_bin := factor(date_bin, levels=c(reference_bin, (1:10)[-reference_bin]))]
date_bins[, reference := FALSE]
date_bins[date_bin == reference_bin, reference := TRUE]

# Get the time of day (10 bins) the sample was taken.
pheno[, date_time_str := paste(attendanceDate, appointmentTime)]
pheno[, date_time_str := gsub(" $", "", date_time_str)]
pheno[appointmentTime != "", date_time := as.POSIXct(date_time_str, format="%d%b%Y %H:%M")]
pheno[, from_midnight := hour(date_time)*60 + minute(date_time)]
pheno[, from_first_draw := from_midnight - min(na.omit(from_midnight)) + 1]
pheno <- pheno[order(from_first_draw)]
pheno[!is.na(from_first_draw), time_bin := ceiling(from_first_draw/(max(from_first_draw)/10))]

# Get information about time bins:
time_bins <- pheno[,.N,by=time_bin]
time_start <- pheno[order(from_midnight), .SD[1], by=time_bin][,.(time_bin, start=from_midnight)]
time_bins <- time_bins[time_start, on=.(time_bin), nomatch=0]
time_end <- time_start[, .(time_bin=time_bin-1, end=start)]
time_last <- pheno[order(from_midnight)][!is.na(time_bin)][.N][, .(time_bin, end=from_midnight)]
time_end <- time_end <- rbind(time_end[-1], time_last)
time_bins <- time_bins[time_end, on=.(time_bin), nomatch=0]
time_bins[, duration := paste(end - start, "minutes")]
time_bins[, start := sprintf("%2d:%02d", start %/% 60, start %% 60)]
time_bins[, end := sprintf("%2d:%02d", end %/% 60, end %% 60)]
time_bins[is.na(time_bin), c("start", "end", "duration") := NA_character_]
time_bins <- time_bins[order(time_bin)]

# Give label to the "NA" level so it's not dropped:
pheno[, time_bin := as.character(time_bin)]
pheno[is.na(time_bin), time_bin := "NA"]
time_bins[, time_bin := as.character(time_bin)]
time_bins[is.na(time_bin), time_bin := "NA"]

# Code factor levels
reference_bin <- time_bins[time_bin != "NA"][which.max(N), time_bin]
pheno[, time_bin := factor(time_bin, levels=c(reference_bin, c(1:10, "NA")[-as.integer(reference_bin)]))]
time_bins[, reference := FALSE]
time_bins[time_bin == reference_bin, reference := TRUE]

# subset phenotype data
pheno <- pheno[, .(IID, age=agePulse, sex=factor(sexPulse), bmi, time_bin, date_bin)]

# To obtain time between blood draw and sample processing, we need an older phenotype
# sheet which has the aliquot processing information
proc_info <- fread("data/INTERVAL/reference_files/aliquot_information.csv")
proc_info[, processedDate := gsub(":.*", "", B_PROCESSED)]
proc_info[, processedTime := gsub("^.*?:", "", B_PROCESSED)]
proc_info[, processedDate := as.IDate(processedDate, format="%d%b%y")]
proc_info[, processedTime := as.ITime(processedTime)]
proc_info[, attendanceDate := as.IDate(attendanceDate, format="%d%b%Y")]
proc_info[, appointmentTime := as.ITime(appointmentTime)]
proc_info[is.na(appointmentTime), appointmentTime := proc_info[!is.na(appointmentTime), median(appointmentTime)]] # median impute missing appointment time (1:55 pm)
proc_info[, minsToProcess := (processedDate - attendanceDate)*24*60 + ((hour(processedTime)*60 + minute(processedTime)) - (hour(appointmentTime)*60 + minute(appointmentTime)))]
proc_info[, daysToProcess := minsToProcess / 60 / 24]
proc_info[omics_idmap, on =.(B2_ALIQUOT_ID=soma4000_gwasQC_bl), IID := Affymetrix_gwasQC_bl]
proc_info <- proc_info[, .(IID, daysToProcess = factor(daysToProcess > 1))]
proc_info <- proc_info[IID %in% unique(raw$IID)]

# Load genetic PCs and extract first 10
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Build table of covariates. 
covar <- gwasqc[omics_idmap, on = .(aliquot_id=soma4000_gwasQC_bl), .(IID=Affymetrix_gwasQC_bl, batch=factor(batch))]
covar <- covar[pcs, on=.(IID)]
covar <- covar[pheno, on = .(IID)]
covar <- covar[proc_info, on = .(IID)]

# Output tables of information
fwrite(date_bins, sep="\t", quote=FALSE, file=sprintf("%s/date_bins.tsv", out_dir))
fwrite(time_bins, sep="\t", quote=FALSE, file=sprintf("%s/time_bins.tsv", out_dir))
fwrite(covar[, .(variable=c("sex", "batch", "days to sample processing"), 
                 alt_group = c("Female", "Batch 2", ">= 1 day"),
                 ref_group = c("Male", "Batch 1", "< 1 day"),
                 alt_N = c(sum(sex == 2), sum(batch == 2), sum(daysToProcess == "TRUE")),
                 ref_N = c(sum(sex == 1), sum(batch == 1), sum(daysToProcess == "FALSE")))],
       sep="\t", quote=FALSE, file=sprintf("%s/factor_covariate_tabulation.tsv", out_dir))

# Summarise the distribution of each raw aptamer measurement
distrib <- raw[, .(
  min = min(value),
  less_5SD = sum(value < (mean(value) - 5*sd(value))),
  l25 = quantile(value, probs=0.25),
  median = median(value),
  mean = mean(value),
  u75 = quantile(value, probs=0.75),
  more_5SD = sum(value > (mean(value) + 5*sd(value))),
  max = max(value)
), by=Seq.ID]

# Build table of associations between raw aptamer measurements and covariates
covar_assocs <- foreach(apt = unique(raw$Seq.ID), .combine=rbind) %dopar% {
  foreach(cvr = names(covar)[-1], .combine=rbind) %do% {
    gc()
    dat <- raw[Seq.ID == apt][covar[, c("IID", cvr), with=FALSE], on=.(IID)]
    l1 <- lm(as.formula(paste0("value ~ ", cvr)), data=dat)
    coef_dt <- as.data.table(coef(summary(l1)), keep.rownames="coefficient")[-1, .(coefficient, Beta=Estimate, Pvalue=`Pr(>|t|)`)]
    ci_dt <- as.data.table(confint(l1), keep.rownames="coefficient")[-1, .(coefficient, L95=`2.5 %`, U95=`97.5 %`)]
    dt <- coef_dt[ci_dt, on = .(coefficient)]
    dt[, Seq.ID := apt]
    dt <- dt[,.(Seq.ID, coefficient, Beta, L95, U95, Pvalue)]
    return(dt)
  }
}

# Cast to wide format
covar_assocs <- dcast(covar_assocs, Seq.ID ~ coefficient, value.var=c("Beta", "L95", "U95", "Pvalue"))

# Add distributions:
covar_assocs <- distrib[covar_assocs, on = .(Seq.ID)]

# Add identifying protein information:
covar_assocs <- prot_info[covar_assocs, on = .(SeqId=Seq.ID)]

# Load associations to provide some sort of ordering:
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- prs_assocs[, .(minFDR = min(FDR)), by=.(Gene, UniProt, Protein=Target)][order(minFDR)]
covar_assocs <- covar_assocs[prs_assocs[,.(Gene, Protein, UniProt)], on = .(Gene, Protein, UniProt)]

# write out
fwrite(covar_assocs, sep="\t", quote=FALSE, file=sprintf("%s/raw_aptamer_covariate_assocs.tsv", out_dir))

# Now, do the same for the post-qc aptamer levels
soma <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")

# Adjust for batch
batch <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")
soma <- soma[batch, on = .(IID), nomatch=0]
soma[, value := lm(value ~ factor(batch))$residuals, by=variable]
soma[, soma_ivt_adj_batch := scale(value)]
soma <- soma[,.(IID, variable, soma_ivt_adj_batch)]

soma <- soma[!prev, on = .(IID)]
soma <- soma[IID %in% unique(hes$IID)]

sinfo <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
soma[sinfo, on = .(variable), variable := SeqId]
soma <- soma[variable %in% prot_info$SeqId]
setnames(soma, c("variable", "soma_ivt_adj_batch"), c("Seq.ID", "value"))

# Summarise the distribution of each adjusted and normaised protein measurement
distrib <- soma[, .(
  min = min(value),
  less_5SD = sum(value < (mean(value) - 5*sd(value))),
  l25 = quantile(value, probs=0.25),
  median = median(value),
  mean = mean(value),
  u75 = quantile(value, probs=0.75),
  more_5SD = sum(value > (mean(value) + 5*sd(value))),
  max = max(value)
), by=Seq.ID]

# Build table of associations between adjusted aptamer measurements and covariates
covar_assocs <- foreach(apt = unique(soma$Seq.ID), .combine=rbind) %dopar% {
  foreach(cvr = names(covar)[-1], .combine=rbind) %do% {
    gc()
    dat <- soma[Seq.ID == apt][covar[, c("IID", cvr), with=FALSE], on=.(IID)]
    l1 <- lm(as.formula(paste0("value ~ ", cvr)), data=dat)
    coef_dt <- as.data.table(coef(summary(l1)), keep.rownames="coefficient")[-1, .(coefficient, Beta=Estimate, Pvalue=`Pr(>|t|)`)]
    ci_dt <- as.data.table(confint(l1), keep.rownames="coefficient")[-1, .(coefficient, L95=`2.5 %`, U95=`97.5 %`)]
    dt <- coef_dt[ci_dt, on = .(coefficient)]
    dt[, Seq.ID := apt]
    dt <- dt[,.(Seq.ID, coefficient, Beta, L95, U95, Pvalue)]
    return(dt)
  }
}

# Cast to wide format
covar_assocs <- dcast(covar_assocs, Seq.ID ~ coefficient, value.var=c("Beta", "L95", "U95", "Pvalue"))

# Add distributions:
covar_assocs <- distrib[covar_assocs, on = .(Seq.ID)]

# Add identifying protein information:
covar_assocs <- prot_info[covar_assocs, on = .(SeqId=Seq.ID)]

# order
covar_assocs <- covar_assocs[prs_assocs[,.(Gene, Protein, UniProt)], on = .(Gene, Protein, UniProt)]

# write out
fwrite(covar_assocs, sep="\t", quote=FALSE, file=sprintf("%s/adjusted_aptamer_covariate_assocs.tsv", out_dir))




