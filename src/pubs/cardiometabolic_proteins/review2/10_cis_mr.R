library(data.table)
library(foreach)
library(doMC)
library(MendelianRandomization)
source("src/utilities/prot_pval.R")

# Set up parallel environment
ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
ncores <- as.integer(ncores)
if (is.na(ncores)) ncores <- 1
registerDoMC(ncores)
setDTthreads(ncores)

# Load IVs
cis_ivs <- fread("analyses/pub/cardiometabolic_proteins/review2/cis_ivs.txt")

# Function to extract relevant coefficients from MR output
extract_coef <- function(mro, egger=FALSE) {
	if (!egger) {
		data.table(mr_estimate = mro@Estimate,
							 mr_se = mro@StdError,
							 mr_L95 = mro@CILower,
							 mr_U95 = mro@CIUpper,
							 mr_pval = mro@Pvalue)
	} else {
		data.table(mr_estimate = c(mro@Estimate, mro@Intercept),
							 mr_se = c(mro@StdError.Int),
							 mr_L95 = c(mro@CILower.Est, mro@CILower.Int),
							 mr_U95 = c(mro@CIUpper.Est, mro@CIUpper.Int),
							 mr_pval = c(mro@Causal.pval, mro@Pleio.pval))
	}
}

# Run MR
tests <- unique(cis_ivs[,.(PRS, GWAS, Target, UniProt, Gene, Aptamer)])
apt_mr <- foreach(test_idx = tests[,.I], .combine=rbind) %dopar% {
  # Get instruments for this aptamer and GWAS
  this_test <- tests[test_idx]
  this_ivs <- cis_ivs[this_test, on = .(PRS, GWAS, Target, UniProt, Gene, Aptamer)]
 
  # Construct MR Input object
  mri <- this_ivs[, mr_input(bx=pQTL.beta, bxse=pQTL.SE, by=gwas.beta, byse=gwas.SE)]

  # If only 1 or 2 instruments can only run IVW
  if (nrow(this_ivs) < 3) {
    ivw <- mr_ivw(mri)
    mro <- rbind(idcol="Method", "IVW"=extract_coef(ivw))
  } else {
		# Run 5 MR methods
		ivw <- mr_ivw(mri)
		simple_median <- mr_median(mri, weighting="simple")
		weighted_median <- mr_median(mri, weighting="weighted")
		weighted_mode <- mr_mbe(mri, weighting="weighted", stderror="simple")
		egger <- mr_egger(mri)

		mro <- rbind(idcol="Method",
			"IVW"=extract_coef(ivw),
			"Median"=extract_coef(simple_median),
			"Weighted Median"=extract_coef(weighted_median),
			"Weighted Mode"=extract_coef(weighted_mode),
			"MR Egger"=extract_coef(egger, TRUE)
		)
		mro[6, Method := "(Intercept)"]
  }

  mro <- cbind(this_test, mro)
  mro
}

# Average across aptamers
prot_mr <- apt_mr[, .(mr_estimate = mean(mr_estimate), mr_se = mean(mr_se), 
                      mr_L95 = mean(mr_L95), mr_U95 = mean(mr_U95), 
                      mr_pval = prot_pvalue(mr_pval, mr_estimate)),
                  by = .(PRS, GWAS, Target, UniProt, Gene, Method)]

# Take median of causal estimates
mr_summary <- prot_mr[Method != "(Intercept)", 
		.(mr_estimate = median(mr_estimate), mr_se = median(mr_se),
			mr_L95 = median(mr_L95), mr_U95 = median(mr_U95),
			mr_pval = prot_pvalue(mr_pval, mr_estimate, median)),
  by = .(PRS, GWAS, Target, UniProt, Gene)]

# FDR adjust
mr_summary[, mr_fdr := p.adjust(mr_pval, method="fdr"), by=.(PRS, GWAS)]

# Row-order
mr_summary <- mr_summary[order(mr_pval)][order(GWAS)][order(PRS)]

# Add pleiotropy p-value
mr_summary[prot_mr[Method == "(Intercept)"], on = .(PRS, GWAS, Target, UniProt, Gene), pleiotropy_pval := i.mr_pval]

# Write out
fwrite(apt_mr, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/mr_all_aptamers.txt")
fwrite(prot_mr, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/mr_all_proteins.txt")
fwrite(mr_summary, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/mr_summary_all_proteins.txt")
