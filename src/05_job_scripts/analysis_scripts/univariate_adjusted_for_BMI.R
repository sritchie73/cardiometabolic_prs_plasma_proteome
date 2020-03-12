#!/usr/bin/env Rscript

# Runs a univariate association scan between a GRS profile
# and an omic platform specified by the input arguments, 
# adjusting for BMI as a covariate
library(data.table)
library(ggplot2)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
grs_name <- args[1]
omic_platform <- args[2]
out_dir <- args[3]

analysis_dir <- file.path(out_dir, grs_name, omic_platform)

grs_profile <- fread(file.path("analyses/GRS_profiles", grs_name, "profile.sscore.gz"), colClasses=c("IID"="character"))
grs_profile <- grs_profile[,.(IID, score)]

trait_measures <- fread(file.path("analyses/processed_traits", omic_platform, "traits.tsv"), colClasses=c("IID"="character"))

# load phenotype data and construct BMI measure - there are some obviously bad
# BMI measures in here due to incorrect weight or height codings, however its 
# not clear how, why, or what cutoffs to use, so I've just chosen clinically
# sensible ones.
pheno <- fread("analyses/processed_traits/phenotypes.tsv", integer64="character", na.strings=c("NA", ""), colClasses=c("IID"="character"))
pheno[, bmi := wt_bl/ht_bl^2]
pheno[wt_bl == 777, bmi := NA_real_] # bad coding
pheno[ht_bl < 1.47, bmi := NA_real_] # clinical cutoff for dwarfism
pheno[ht_bl > 2.1, bmi := NA_real_] # clinical cutoff for gigantism
pheno[wt_bl < 50 | wt_bl > 160, bmi := NA_real_] # NHS restrictions for weight

# Load PCs to use as covariates
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt", colClasses=c("ID"="character"))
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

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

covariates <- merge(covariates, pheno[,.(IID, bmi)], by = "IID", all = FALSE)

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

