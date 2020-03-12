#!/usr/bin/env Rscript

# Runs a univariate association scan, adjusting for date as a covariate.
# Since seasonality effects are not necessarily linear with time, we treat this
# as a catagorical variable. We split the attendanceDate variable into ten
# bins of equal duration, restricting to samples in the platform of interest.
# In the full dataset this results in ten bins of 73 - 74 days:
#
#     date_bin    N      start        end duration
#  1:        1 4105 2012-06-11 2012-08-23       73
#  2:        2 9641 2012-08-23 2012-11-05       74
#  3:        3 7670 2012-11-05 2013-01-17       73
#  4:        4 6205 2013-01-17 2013-04-01       74
#  5:        5 5421 2013-04-01 2013-06-13       73
#  6:        6 3976 2013-06-13 2013-08-26       74
#  7:        7 4486 2013-08-26 2013-11-06       72
#  8:        8 3697 2013-11-06 2014-01-19       74
#  9:        9 2374 2014-01-19 2014-04-02       73
# 10:       10 1097 2014-04-02 2014-06-14       73
#
# We treat each platform separately for the purposes of date binning because their
# samples may not cover the full range of dates. For example the samples with 
# somalogic protein measurements fall between 2012-06-15 to 2013-10-23, so we'd 
# only get bins 1-7 in that dataset. The exact binning for the platform is output
# along with the adjusted associations.
#
# Note that each bin contains a different number of samples. However, splitting 
# the data into equal time duration bins likely better reflects any seasonality 
# effects on protein levels (or metabolites etc.). In the protein data, the 
# median protein level in each of these bins more closely matched the loess fit 
# across all data points with date of blood draw as the dependent variable when 
# compared to binning the data into ten equal sized sample bins. 
# See src/sanity_checks/circadian_seasonality_effects.R and
# analyses/sanity_checks/seasonality_effects.*.png
library(data.table)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
grs_name <- args[1]
omic_platform <- args[2]
out_dir <- args[3]

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

# Load phenotype data and construct date bins. Subset to omic platform samples to do so.
pheno <- fread("analyses/processed_traits/phenotypes.tsv", integer64="character", na.strings=c("NA", ""))
pheno[, date := as.POSIXct(attendanceDate, format="%d%b%Y")]
pheno <- pheno[,.(IID, date)]
pheno <- pheno[IID %in% trait_measures$IID]
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

# Add to covariates
if (exists("covariates")) {
  covariates <- merge(covariates, pheno[,.(IID, date_bin)], by="IID")
} else {
  covariates <- pheno[,.(IID, date_bin)]
} 

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
fwrite(date_bins, file=file.path(analysis_dir, "date_bins.tsv"), sep="\t", quote=FALSE)

