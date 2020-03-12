# Retest associations
library(data.table)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)
pairs <- fread(args[1], header=FALSE)
setnames(pairs, c("grs", "variable"))

proteins <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv", colClasses=c("IID"="character"))

# Load PCs to use as covariates
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt", colClasses=c("ID"="character"))
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Retest each pair
for (ii in seq_len(nrow(pairs))) {
  grs <- pairs[ii, grs]
  aptamer <- pairs[ii, variable]
  
  # Load protein specific and GRS specific covariates
  platform_covar_file <- file.path("analyses/processed_traits/somalogic_proteins/covariates.tsv")
  grs_covar_file <- file.path("analyses/GRS_profiles/", grs, "covariates.tsv")
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

	# Build the model formula:
  model <- paste("value ~ score +", paste(names(covariates)[-1], collapse=" + "))
  model <- gsub("gender", "factor(gender)", model)
  model <- gsub("batch", "factor(batch)", model)
  model <- gsub("samplegroup", "factor(samplegroup)", model)

  # Load GRS levels re-calculated after removing this aptamers pQTLs:
  levels <- fread(sprintf("analyses/grs_pqtl_removed/%s/%s/profile.sscore.gz", grs, aptamer), colClasses=c("IID"="character"))

  # Create combined dataset for model testing
  dat <- merge(levels[, .(IID, score)], proteins[variable == aptamer], by="IID")
  dat <- merge(dat, covariates, by="IID")

  # Fit model and extract coefficients
  l1 <- lm(as.formula(model), data=dat)
  l1_coef <- coef(summary(l1))
  ci95 <- confint(l1)

  # Create table of results and save to directory
  res <- data.table(grs=grs, trait=aptamer, beta=l1_coef[2,1], 
                    l95=ci95[2,1], u95=ci95[2,2], pval=l1_coef[2,4])
  fwrite(res, sep="\t", quote=FALSE, 
    file=sprintf("analyses/grs_pqtl_removed/%s/%s/associations.tsv", grs, aptamer))
}


