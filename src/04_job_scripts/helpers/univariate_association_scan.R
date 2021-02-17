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
  data.table(trait=trait, beta=l1_coef[2,1], l95=ci95[2,1], u95=ci95[2,2], pval=l1_coef[2,4],
             tstat = l1_coef[2,3], fstat = l1_summary$fstatistic[1], df=l1_summary$fstatistic[3], r2=l1_summary$adj.r.squared)
}

# Run stratified analysis, if requested
strata_file <- sprintf("analyses/GRS_profiles/%s/stratify.txt", grs_name)
if (file.exists(strata_file)) {
  strata <- fread(strata_file)
  strata <- melt(strata, id.vars="IID", variable.name="group_name", value.name="group") 

  # Stratify by each group variable separately 
  strata_res <- foreach(gp = strata[,unique(group_name)], .combine=rbind) %do% {
    foreach(gpv = strata[group_name == gp, unique(group)], .combine=rbind) %do% { # Each group value, e.g. Sex == 1, Sex == 2
      foreach(trait = unique(trait_measures$variable), .combine=rbind) %do% {

        # Remove group variable from model, e.g. if its a covariate
        if (grepl(sprintf(" %s ", gpv), model))
          model <- gsub(sprintf("+ %s ", gpv), "", model)
        if (grepl(sprintf("factor(%s)", gpv), model))
          model <- gsub(sprintf("+ factor(%s)", gpv), "", model)

        # build data for linear model
        dat <- merge(trait_measures[variable == trait], covariates, by = "IID")
        dat <- merge(grs_profile, dat, by = "IID")
        dat <- dat[IID %in% strata[group_name == gp & group == gpv, IID]]

        # Fit model and extract results
        l1 <- lm(as.formula(model), data=dat)
        l1_summary <- summary(l1)
        l1_coef <- coef(l1_summary)
        ci95 <- confint(l1)
        data.table(group_name = gp, group = gpv, trait=trait, beta=l1_coef[2,1], l95=ci95[2,1], u95=ci95[2,2], pval=l1_coef[2,4],
                   tstat = l1_coef[2,3], fstat = l1_summary$fstatistic[1], df=l1_summary$fstatistic[3], r2=l1_summary$adj.r.squared)
      }
    } 
  }
  res[, c("group_name", "group") := .("all", NA)]
  res <- rbind(res, strata_res, use.names=TRUE)

  res <- res[order(abs(beta))][order(pval)][order(group)][order(group_name)]
  res[, fdr := p.adjust(pval, method="fdr"), by=.(group_name, group)]
} else {
  res <- res[order(abs(beta))][order(pval)] 
  res[, fdr := p.adjust(pval, method="fdr")]
}

# Write out table of results
fwrite(res, file=file.path(analysis_dir, "associations.tsv"), sep="\t", quote=FALSE)
