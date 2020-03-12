#!/usr/bin/env Rscript

# Runs a univariate association scan between a GRS profile
# and an omic platform specified by the input arguments, 
# adjusting for any significant QTLs as covariates
library(data.table)
library(ggplot2)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
grs_name <- args[1]
omic_platform <- args[2]
out_dir <- args[3]

analysis_dir <- file.path(out_dir, grs_name, omic_platform)
omic_dir <- file.path("analyses/processed_traits", omic_platform)

if (!file.exists(file.path(omic_dir, "qtl_list.txt"))) {
  # No QTLs for this platform, prune empty directory and exit
  system(sprintf("rmdir %s", analysis_dir))
  quit(save="no", status=0)
}

grs_profile <- fread(file.path("analyses/GRS_profiles", grs_name, "profile.sscore.gz"), colClasses=c("IID"="character"))
grs_profile <- grs_profile[,.(IID, score)]
trait_measures <- fread(file.path("analyses/processed_traits", omic_platform, "traits.tsv"), colClasses=c("IID"="character"))

# Load list of QTL-trait associations
qtl_assocs <- fread(sprintf("%s/qtl_list.txt", omic_dir))

# Load in genotype data
qtl_geno <- fread(sprintf("%s/qtl_geno_prob.tsv", omic_dir), colClasses=c("IID"="character"))

# We can't have variant ids of the form chr:pos here, because : is a 
# special reserved symbol in model formulae
colnames(qtl_geno) <- gsub(":", "_", colnames(qtl_geno))
qtl_assocs[, variant := gsub(":", "_", variant)]

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

# Determine whether the trait ~ GRS association is attenuated when adjusting for QTLs
base_model <- paste("value ~ score +", paste(names(covariates)[-1], collapse=" + "))
base_model <- gsub("gender", "factor(gender)", base_model)
base_model <- gsub("batch", "factor(batch)", base_model)
base_model <- gsub("samplegroup", "factor(samplegroup)", base_model)
base_model <- gsub("time_code", "factor(time_code)", base_model)

res <- foreach(trait = unique(trait_measures$variable), .combine=rbind) %do% {
  X <- merge(trait_measures[variable == trait], covariates, by = "IID")
  dat <- merge(grs_profile, X, by = "IID")
  
  # Add QTLs for that trait if they exist
  model <- base_model
  if (trait %in% qtl_assocs$variable) {
    qtl_vars <- qtl_assocs[variable == trait, variant]
    qtl_vars <- intersect(qtl_vars, colnames(qtl_geno)) # filter for QTLs that didnt pass QC in INTERVAL.
    if (length(qtl_vars) > 0) {
			dat <- merge(dat, qtl_geno[, c("IID", qtl_vars), with=FALSE], by="IID")
			model <- paste(base_model, "+", paste(qtl_vars, collapse=" + "))
    } else {
      stop(trait)
    }
  }

  l1 <- lm(as.formula(model), data=dat)

  # Extract model results
  l1_coef <- coef(summary(l1))
  ci95 <- confint(l1)
 
	data.table(trait=trait, beta=l1_coef[2,1], 
						 l95=ci95[2,1], u95=ci95[2,2], 
             pval=l1_coef[2,4])
}
res[, fdr := p.adjust(pval, method="fdr")]
fwrite(res, file=file.path(analysis_dir, "associations.tsv"), sep="\t", quote=FALSE)

# Determine the contributions of the GRS compared to the QTLs for the trait in question
res <- foreach(trait = unique(trait_measures$variable), .combine=rbind) %do% {
  X <- merge(trait_measures[variable == trait], covariates, by = "IID")
  dat <- merge(grs_profile, X, by = "IID")
  
  # Add QTLs for that trait if they exist
  model <- base_model
  if (trait %in% qtl_assocs$variable) {
    qtl_vars <- qtl_assocs[variable == trait, variant]
    qtl_vars <- intersect(qtl_vars, colnames(qtl_geno)) # filter for QTLs that didnt pass QC in INTERVAL.
    if (length(qtl_vars) > 0) {
			dat <- merge(dat, qtl_geno[, c("IID", qtl_vars), with=FALSE], by="IID")
			model <- paste(base_model, "+", paste(qtl_vars, collapse=" + "))
    } else {
      stop(trait)
    }
  }

  l1 <- lm(as.formula(model), data=dat)

  # Extract model results
  l1_coef <- coef(summary(l1))
  ci95 <- confint(l1)
 
  # Make sure we extract the qtl coefficients as well
  if (trait %in% qtl_assocs$variable) {
    coef_names <- c("score", qtl_vars)
    coef_dt <- as.data.table(l1_coef, keep.rownames="coef")
    data.table(trait=trait, variable=c("grs", coef_names[-1]),
               beta=coef_dt[coef_names, on = .(coef), Estimate],
               l95=ci95[coef_names,1], u95=ci95[coef_names,2],
               pval=coef_dt[coef_names, on = .(coef), `Pr(>|t|)`])
  } else {
		data.table(trait=trait, variable="grs", 
							 beta=l1_coef[2,1], 
							 l95=ci95[2,1], 
							 u95=ci95[2,2], 
							 pval=l1_coef[2,4])
  }
}
fwrite(res, file=file.path(analysis_dir, "qtl_grs_associations.tsv"), sep="\t", quote=FALSE)
