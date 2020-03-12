#!/usr/bin/env Rscript

# Runs a univariate association scan, adjusting for time of day as a covariate.
# Since Circadian effects are not necessarily linear with time, we treat this
# as a catagorical variable. We split the appointmentTime variable into 
# ten bins of equal duration from the time of first blood draw on any day 
# to the time of last blood draw of any day. This is done within each platform
# separately, in case (for whatever reason) the fall range of times is not
# covered in platform-specific sample subset and we lose some bins. In the full
# dataset the 10 bins and durations look as follows:
#
#     time_bin    N start   end   duration
#  1:        1 1379  7:35  8:55 80 minutes
#  2:        2 3744  8:55 10:10 75 minutes
#  3:        3 4751 10:10 11:25 75 minutes
#  4:        4 5367 11:25 12:45 80 minutes
#  5:        5 5657 12:45 14:00 75 minutes
#  6:        6 6136 14:00 15:15 75 minutes
#  7:        7 4780 15:15 16:35 80 minutes
#  8:        8 3884 16:35 17:50 75 minutes
#  9:        9 3668 17:50 19:05 75 minutes
# 10:       10  985 19:05 20:20 75 minutes
# 11:       NA 8321  <NA>  <NA>       <NA>
#
# In the somalogic protein sample subset the range of possible times is from 8am
# to 8:15 pm. The specific bins for each platform are output along with the 
# adjusted associations.
#
# Note that each bin contains a different number of samples. However, splitting 
# the data into equal time duration bins likely better reflects any time of day
# effects on protein levels (or metabolites etc.). In the protein data, the 
# median protein level in each of these bins more closely matched the loess fit 
# across all data points with time of blood draw (hour and minute) as the 
# dependent variable when compared to binning the data into ten equal sized 
# sample bins. See src/sanity_checks/circadian_seasonality_effects.R and
# analyses/sanity_checks/circadian_effects.*.png
#
# Note also a large proportion of samples are missing appointmentTime information
# (N=480/3175 for the protein data), these are excluded in the sensitivity analysis.
library(data.table)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
grs_name <- args[1]
omic_platform <- args[2]
out_dir <- args[3]

if (omic_platform == "metabolon_metabolomics") {
  quit(save="no", status=0) # already adjusted for time of day, nothing to do here
}

analysis_dir <- file.path(out_dir, grs_name, omic_platform)

# Load GRS levels
grs_profile <- fread(file.path("analyses/GRS_profiles", grs_name, "profile.sscore.gz"), colClasses=c("IID"="character"))
grs_profile <- grs_profile[,.(IID, score)]

# Load omic measurements
trait_measures <- fread(file.path("analyses/processed_traits", omic_platform, "traits.tsv"), colClasses=c("IID"="character"))

# Load PCs to use as covariates
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt", colClasses=c("ID"="character"))
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Load covariates for both the omics platform and GRS
platform_covar_file <- file.path("analyses/processed_traits", omic_platform, "covariates.tsv")
grs_covar_file <- file.path("analyses/GRS_profiles/", grs_name, "covariates.tsv")
if (file.exists(platform_covar_file) && !file.exists(grs_covar_file)) {
  covariates <- fread(platform_covar_file, colClasses=c("IID"="character"))
} else if (file.exists(grs_covar_file) && !file.exists(platform_covar_file)) {
  covariates <- fread(grs_covar_file, colClasses=c("IID"="character"))
} else if (file.exists(grs_covar_file) && file.exists(platform_covar_file)) {
  platform_covariates <- fread(platform_covar_file, colClasses=c("IID"="character"))
  grs_covariates <- fread(grs_covar_file, colClasses=c("IID"="character"))
  covariates <- merge(platform_covariates, grs_covariates, by = "IID")
}

if (exists("covariates")) {
  covariates <- merge(covariates, pcs, by="IID")
} else {
    covariates <- pcs
}

# Load phenotype data and construct time series bins:
pheno <- fread("analyses/processed_traits/phenotypes.tsv", integer64="character", na.strings=c("NA", ""))
pheno[, date_time_str := paste(attendanceDate, appointmentTime)]
pheno[, date_time_str := gsub(" $", "", date_time_str)]
pheno[appointmentTime != "", date_time := as.POSIXct(date_time_str, format="%d%b%Y %H:%M")]
pheno <- pheno[,.(IID, date_time)]

pheno <- pheno[IID %in% trait_measures$IID]
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

# Code factor levels
reference_bin <- time_bins[!is.na(time_bin)][which.max(N), time_bin]
pheno[, time_bin := factor(time_bin, levels=c(reference_bin, (1:10)[-reference_bin]))]
time_bins[, reference := FALSE]
time_bins[time_bin == reference_bin, reference := TRUE]

# Add to covariates
covariates <- merge(covariates, pheno[,.(IID, time_bin)], by="IID")

# Build the model formula:
model <- paste("value ~ score +", paste(names(covariates)[-1], collapse=" + "))
model <- gsub("gender", "factor(gender)", model)
model <- gsub("batch", "factor(batch)", model)
model <- gsub("samplegroup", "factor(samplegroup)", model)
model <- gsub("time_code", "factor(time_code)", model)

# Run for each trait
res <- foreach(trait = unique(trait_measures$variable), .combine=rbind) %do% {
  X <- merge(trait_measures[variable == trait], covariates, by = "IID")
  dat <- merge(grs_profile, X, by = "IID")
  l1 <- lm(as.formula(model), data=dat)

  # Extract model results
  l1_coef <- coef(summary(l1))
  ci95 <- confint(l1)

  data.table(trait=trait, beta=l1_coef[2,1], l95=ci95[2,1], u95=ci95[2,2], pval=l1_coef[2,4])
}
res <- res[order(abs(beta))][order(pval)] 

# FDR correction
res[, fdr := p.adjust(pval, method="fdr")]

# Write out table of results
fwrite(res, file=file.path(analysis_dir, "associations.tsv"), sep="\t", quote=FALSE)

# Write out bin information
fwrite(time_bins, file=file.path(analysis_dir, "time_bins.tsv"), sep="\t", quote=FALSE)

