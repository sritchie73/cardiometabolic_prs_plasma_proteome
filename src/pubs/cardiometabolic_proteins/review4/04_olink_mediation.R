library(data.table)
library(foreach)
library(doMC)
library(survival)
library(lubridate)
library(medflex)

source("src/utilities/prot_pval.R")

# Set up parallel environment
if (!exists("ncores")) {
  ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
  ncores <- as.integer(ncores)
  if(is.na(ncores)) ncores <- 1
}
registerDoMC(ncores)
setDTthreads(ncores)

# Load HES data
source("src/pubs/cardiometabolic_proteins/review2/load_HES.R")

# Filter to primary+secondary endpoints
hes <- hes[code_type == "primary+secondary"]

# Load list of cardiometabolic events for prevalent event exclusion
source("src/pubs/cardiometabolic_proteins/review2/cardiometabolic_events.R")

# Exclude people with prevalent events
prev_24m <- unique(hes[prevalent_24m == 1 & phenotype %in% cardiometabolic, .(IID)])

# Load PRS
prs <- rbind(idcol="PRS",
  CAD_PRS = fread("analyses/GRS_profiles/CAD_metaGRS/profile.sscore.gz"),
  IS_PRS = fread("analyses/GRS_profiles/Stroke_metaGRS/profile.sscore.gz"),
  T2D_PRS = fread("analyses/GRS_profiles/T2D_2018/profile.sscore.gz"),
  AF_PRS = fread("analyses/GRS_profiles/Afib_2018/profile.sscore.gz"),
  CKD_PRS = fread("analyses/GRS_profiles/CKD_2019/profile.sscore.gz")
)

# Load PCs
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Adjust PRS for PCs
prs <- prs[pcs, on = .(IID), nomatch=0]
prs[, prs_adj_pcs := lm(scale(score_sum) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +
                                           PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=PRS]
prs <- prs[, .(IID, PRS, prs_adj_pcs)]

# Load phenotype data
pheno <- fread("analyses/processed_traits/phenotypes.tsv")
pheno <- pheno[!is.na(IID)]

# Compute age at follow-up
pheno <- pheno[,.(IID, sex = sexPulse, age_bl = agePulse, ht_bl, wt_bl, attendanceDate, attendanceDate_24m)]
pheno[, attendanceDate := as.IDate(attendanceDate, format="%d%B%Y")]
pheno[, attendanceDate_24m := as.IDate(attendanceDate_24m, format="%d%B%Y")]
pheno[, follow_24m := time_length(as.Date(attendanceDate_24m) - as.Date(attendanceDate), unit="year")]
pheno[, age_24m := age_bl + round(follow_24m, digits=1)]

# Filter columns
pheno <- pheno[,.(IID, sex, age_24m)]

# Flag prevalent cases
pheno[, prev_24m := ifelse(is.na(age_24m), NA, FALSE)]
pheno[prev_24m, on = .(IID), prev_24m := TRUE]

# Load SomaLogic aptamer levels (so we know who to exclude)
soma <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")

# Load in olink qPCR data
olink_qpcr <- fread("analyses/processed_traits/olink_proteins/traits.tsv")

# Remove people with prevalent cardiometabolic disease at the olink timepoint
olink_qpcr <- olink_qpcr[!prev_24m, on = .(IID)]

# Remove samples with somalogic proteomics
olink_qpcr <- olink_qpcr[!(IID %in% unique(soma$IID))]

# Load in PRS to protein associations and filter to those with FDR < 0.05
soma_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
soma_assocs <- soma_assocs[FDR < 0.05]

# Load olink information sheet and filter to proteins also measured on the somalogic platform
olink_qpcr_info <- fread("analyses/processed_traits/olink_proteins/trait_info.tsv")
olink_qpcr_info <- olink_qpcr_info[panel != "neu"] # agreement not in place to use Neurological panel.
common_prot <- intersect(unique(soma_assocs$UniProt), olink_qpcr_info$UniProt)

# Get common tests
common_tests <- unique(soma_assocs[UniProt %in% common_prot, .(PRS, UniProt, Gene)])

# Filter endpoint data
hes <- hes[phenotype %in% c("Myocardial infarction", "Diabetes")]
common_tests[PRS == "CAD_PRS", endpoint := "Myocardial infarction"]
common_tests[PRS == "T2D_PRS", endpoint := "Diabetes"]

# Get HRs for incident disease:
hrs <- foreach(tidx = common_tests[,.I], .combine=rbind) %do% {
	this_prs <- common_tests[tidx, PRS]
	this_up <- common_tests[tidx, UniProt]
	this_gene <- common_tests[tidx, Gene]
  this_var <- olink_qpcr_info[UniProt == this_up, variable]
  this_end <- common_tests[tidx, endpoint]

  dat <- olink_qpcr[variable == this_var]
  dat <- dat[pheno, on = .(IID), nomatch=0]
  dat <- dat[hes[phenotype == this_end], on = .(IID), nomatch=0]

	c1 <- coxph(Surv(followUp_24m, event) ~ scale(value) + age_24m + factor(sex), data=dat)
  cf <- coef(summary(c1))
  ci <- confint(c1)

	data.table(samples = dat[,.N], endpoint = this_end, cases=dat[,sum(event)], 
             med_onset_follow = dat[event == 1, median(followUp_24m)], med_onset_age = dat[event == 1, median(age_24m)],
             med_follow = dat[,median(followUp_24m)], max_follow = dat[,max(followUp_24m)],
             Gene = this_gene, logHR = cf[1,1], logHR.SE = cf[1,3], HR = cf[1,2], HR.L95 = exp(ci[1,1]),
					   HR.U95 = exp(ci[1,2]), HR.P = cf[1,5])
}

# Mediation analysis
med <- foreach(tidx = common_tests[,.I], .combine=rbind) %do% {
  this_prs <- common_tests[tidx, PRS]
  this_up <- common_tests[tidx, UniProt]
  this_gene <- common_tests[tidx, Gene]
  this_var <- olink_qpcr_info[UniProt == this_up, variable]
  this_end <- common_tests[tidx, endpoint]

  dat <- olink_qpcr[variable == this_var]
  dat <- dat[pheno, on = .(IID), nomatch=0]
  dat <- dat[hes[phenotype == this_end], on = .(IID), nomatch=0]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]
  dat[, value := scale(value)]
  dat[, prs_adj_pcs := scale(prs_adj_pcs)]
  dat[, age_24m := scale(age_24m)]
  dat[, sex := factor(sex)]
  dat[, event := factor(event)]

  extData <- neImpute(event ~ prs_adj_pcs + value + age_24m + sex, family = binomial("logit"), data=dat)
  nMod <- neModel(event ~ prs_adj_pcs0 + prs_adj_pcs1 + age_24m + sex, family = binomial("logit"), expData=extData, se="robust")
  cf <- coef(summary(neEffdecomp(nMod)))
  ci <- confint(neEffdecomp(nMod))

  data.table(PRS = this_prs, Endpoint = this_end, Gene = this_gene, Mediation = rownames(cf), 
             logOR = cf[,1], SE = cf[, 2], L95 = ci[,1], U95 = ci[,2], P = cf[,4])
}

# split into effect types
nde <- med[Mediation == "natural direct effect"]
nie <- med[Mediation == "natural indirect effect"]
te <- med[Mediation == "total effect"]

# Compute % mediated and 95% CIs
nde[te, on = .(PRS, Gene, Endpoint), PTE := logOR / i.logOR]
nde[te, on = .(PRS, Gene, Endpoint), PTE.L95 := L95 / i.logOR]
nde[te, on = .(PRS, Gene, Endpoint), PTE.U95 := U95 / i.logOR]
nde <- nde[, .(PRS, Gene, Endpoint, Mediation, logOR, SE, L95, U95, PTE, PTE.L95, PTE.U95, P)]

nie[te, on = .(PRS, Gene, Endpoint), PTE := logOR / i.logOR]
nie[te, on = .(PRS, Gene, Endpoint), PTE.L95 := L95 / i.logOR]
nie[te, on = .(PRS, Gene, Endpoint), PTE.U95 := U95 / i.logOR]
nie <- nie[, .(PRS, Gene, Endpoint, Mediation, logOR, SE, L95, U95, PTE, PTE.L95, PTE.U95, P)]

te <- te[, .(PRS, Gene, Endpoint, Mediation, logOR, SE, L95, U95, PTE=1, PTE.L95=1, PTE.U95=1, P)]

# Combine
med <- rbind(nde, nie, te)

# Compute PGS OR
common_tests <- unique(common_tests[,.(PRS, endpoint)])
ors <- foreach(tidx = common_tests[,.I], .combine=rbind) %do% {
  this_prs <- common_tests[tidx, PRS]
  this_end <- common_tests[tidx, endpoint]

  dat <- unique(olink_qpcr[variable %in% olink_qpcr_info[UniProt %in% common_prot, variable], .(IID)])
  dat <- pheno[dat, on = .(IID), nomatch=0]
  dat <- dat[hes[phenotype == this_end], on = .(IID), nomatch=0]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]
  dat[, prs_adj_pcs := scale(prs_adj_pcs)]
  dat[, age_24m := scale(age_24m)]
  dat[, sex := factor(sex)]
  dat[, event := factor(event)]

  c1 <- glm(event ~ prs_adj_pcs + age_24m + factor(sex), family="binomial", data=dat)
  cf <- coef(summary(c1))
  ci <- confint(c1)

  data.table(samples = dat[,.N], endpoint = this_end, cases=dat[event == 1, .N], PRS = this_prs,
             logOR = cf[2,1], logOR.SE = cf[2,2], OR = exp(cf[2,1]), OR.L95 = exp(ci[2,1]),
             OR.U95 = exp(ci[2,2]), OR.P = cf[2,4])
}

# Write out
fwrite(hrs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/olink_hr_replication.txt")
fwrite(med, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/olink_med_replication.txt")
fwrite(ors, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/olink_pgs_disease_assocs.txt")




