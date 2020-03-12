#!/usr/bin/env Rscript

# Runs a univariate association scan between a GRS profile
# and an omic platform specified by the input arguments
library(data.table)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
grs_name <- args[1]
omic_platform <- args[2]

analysis_dir <- file.path("analyses/univariate_associations", grs_name, omic_platform)

grs_profile <- fread(file.path("analyses/GRS_profiles", grs_name, "profile.sscore.gz"))
grs_profile <- grs_profile[,.(IID, score)]

# Load PCs to use as covariates
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

trait_measures <- fread(file.path("analyses/processed_traits", omic_platform, "traits.tsv"))
platform_covar_file <- file.path("analyses/processed_traits", omic_platform, "covariates.tsv")
grs_covar_file <- file.path("analyses/GRS_profiles/", grs_name, "covariates.tsv")
if (file.exists(platform_covar_file) && !file.exists(grs_covar_file)) {
  covariates <- fread(platform_covar_file)
} else if (file.exists(grs_covar_file) && !file.exists(platform_covar_file)) {
  covariates <- fread(grs_covar_file)
} else if (file.exists(grs_covar_file) && file.exists(platform_covar_file)) {
  platform_covariates <- fread(platform_covar_file)
  grs_covariates <- fread(grs_covar_file)
  covariates <- merge(platform_covariates, grs_covariates, by = "IID")
}

if (exists("covariates")) {
  covariates <- merge(covariates, pcs, by="IID")
} else {
  covariates <- pcs
}

# Build the model formula:
model <- paste("value ~ score +", paste(names(covariates)[-1], collapse=" + "))
model <- gsub("gender", "factor(gender)", model)
model <- gsub("batch", "factor(batch)", model)
model <- gsub("samplegroup", "factor(samplegroup)", model)
model <- gsub("time_code", "factor(time_code)", model)

# Run for each trait
res <- foreach(trait = unique(trait_measures$variable), .combine=rbind) %do% {
  # subset data and run linear model
  X <- merge(trait_measures[variable == trait], covariates, by = "IID")
  dat <- merge(grs_profile, X, by = "IID")
  
  # Fit model and extract results
  l1 <- lm(as.formula(model), data=dat)
  l1_summary <- summary(l1)
  l1_coef <- coef(l1_summary)
  ci95 <- confint(l1)
  data.table(trait=trait, beta=l1_coef[2,1], l95=ci95[2,1], u95=ci95[2,2], pval=l1_coef[2,4])
}

# Write out table of results
fwrite(res, file=file.path(analysis_dir, "associations.tsv"), sep="\t", quote=FALSE)
