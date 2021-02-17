library(data.table)
library(foreach)
library(doMC)
library(survival)
library(medflex)
library(mma)
library(ggplot2)
library(scales)
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
prev <- unique(hes[prevalent == 1 & phenotype %in% cardiometabolic, .(IID)])
hes <- hes[!prev, on = .(IID)]

# Load PRS
prs <- rbind(idcol="PRS",
  CAD_PRS = fread("analyses/GRS_profiles/CAD_metaGRS/profile.sscore.gz"),
  IS_PRS = fread("analyses/GRS_profiles/Stroke_metaGRS/profile.sscore.gz"),
  T2D_PRS = fread("analyses/GRS_profiles/T2D_2018/profile.sscore.gz"),
  CKD_PRS = fread("analyses/GRS_profiles/CKD_2019/profile.sscore.gz"),
  AF_PRS = fread("analyses/GRS_profiles/Afib_2018/profile.sscore.gz")
)

# Load PCs
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Adjust PRS for PCs
prs <- prs[pcs, on = .(IID), nomatch=0]
prs[, prs_adj_pcs := lm(scale(score_sum) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +
                                           PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=PRS]
prs <- prs[, .(IID, PRS, prs_adj_pcs)]

# Load phenotype data - use the HES data because it has all samples mapped
pheno_file <- list.files("data/INTERVAL/HES", pattern="INTERVALdata_[^(p|P)3].*", full.names=TRUE)
pheno <- fread(pheno_file)

idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
idmap <- fread(idmap_file)

pheno[idmap, on = .(identifier), IID := Affymetrix_gwasQC_bl]
pheno <- pheno[!is.na(IID), .(IID, age=agePulse, sex=sexPulse, weight=wt_bl, height=ht_bl)]

pheno <- pheno[!prev, on = .(IID)]

# Compute BMI
pheno[weight == 777, weight := NA] # bad coding
pheno[, bmi := weight/height^2]
pheno[height < 1.47, bmi := NA_real_] # clinical cutoff for dwarfism
pheno[height > 2.1, bmi := NA_real_] # clinical cutoff for gigantism
pheno[weight < 50 | weight > 160, bmi := NA_real_] # NHS restrictions for weight

# Load SomaLogic aptamer levels
soma <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")

# Adjust for batch
batch <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")
soma <- soma[batch, on = .(IID), nomatch=0]
soma[, value := lm(value ~ factor(batch))$residuals, by=variable]
soma[, soma_ivt_adj_batch := scale(value)]
soma <- soma[,.(IID, variable, soma_ivt_adj_batch)]

# Filter prs table to samples with somalogic measures
prs <- prs[IID %in% unique(soma$IID)]
prs[, prs_adj_pcs := scale(prs_adj_pcs)] # standardise to SD = 1, mean = 0

# Drop prevalent cases
prs <- prs[!prev, on = .(IID)]
soma <- soma[!prev, on = .(IID)]

# And withdrawn samples
idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
idmap <- fread(idmap_file)

prs <- prs[IID %in% idmap$Affymetrix_gwasQC_bl]
soma <- soma[IID %in% idmap$Affymetrix_gwasQC_bl]
pheno <- pheno[IID %in% idmap$Affymetrix_gwasQC_bl]

# Load Aptamer information
sinfo <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
sinfo <- sinfo[, .(variable, SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, Gene=Gene.Name)]

############################ --- Finished loading data

dat <- NULL # variable needs to be globally defined for neImpute to find it

############################ Test Mediation of T2D PGS to protein associations 

# Load set of potentially BMI mediated T2D PGS to protein associations
bmi_adj <- fread("analyses/pub/cardiometabolic_proteins/review2/bmi_attenuated_assocs.txt")
tests <- bmi_adj[,.(PRS, Target, UniProt, Gene)]
tests <- tests[sinfo, on = .(Target, UniProt, Gene), nomatch=0]

mediation <- foreach(test_idx = tests[,.I], .combine=rbind, .export=c("soma", "prs", "pheno", "tests")) %do% {
  this_prs = tests[test_idx, PRS]
  this_target = tests[test_idx, Target]
  this_uniprot = tests[test_idx, UniProt]
  this_gene = tests[test_idx, Gene]
  this_aptvar = tests[test_idx, variable]
  this_aptamer = tests[test_idx, Aptamer]

  dat <- soma[variable == this_aptvar]
  dat <- dat[pheno, on = .(IID), nomatch=0]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]
  dat[, prs_adj_pcs := scale(prs_adj_pcs)]
  dat[, soma_ivt_adj_batch := scale(soma_ivt_adj_batch)]
  dat[, bmi := scale(log(bmi))]
  dat <- dat[!is.na(bmi)]

	extData <- neImpute(soma_ivt_adj_batch ~ prs_adj_pcs  + bmi, data=dat)
	nMod <- neModel(soma_ivt_adj_batch ~ prs_adj_pcs0 + prs_adj_pcs1, expData=extData, se="robust")
	cf <- coef(summary(neEffdecomp(nMod)))
	ci <- confint(neEffdecomp(nMod))

  data.table(PRS = this_prs, Target = this_target, 
						 UniProt = this_uniprot, Gene = this_gene, Aptamer = this_aptamer,
						 Mediation = rownames(cf), Beta = cf[,1], SE = cf[, 2], L95 = ci[,1],
						 U95 = ci[,2], P = cf[,4], PTE = cf[,1] / cf["total effect", 1]) 
}

# Summarise to protein level
mediation_prot <- mediation[, .(Beta = mean(Beta), SE = mean(SE), L95 = mean(L95), U95 = mean(U95), P = prot_pvalue(P, Beta), PTE = mean(PTE)),
                             by = .(PRS, Target, UniProt, Gene, Mediation)]

# Collate
mediation <- merge(mediation_prot, mediation, by=c("PRS", "Target", "UniProt", "Gene", "Mediation"), suffixes=c("", ".aptamer"))
fwrite(mediation, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/bmi_mediation.txt")

############################ Test Mediation of PGS to disease associations by protein levels

# Load significant PRS to protein associations
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[FDR < 0.05,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])

# Map PGS to endpoints
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  endpoint = c("Atrial fibrillation", "Myocardial infarction", "End stage renal disease", "Ischaemic stroke", "Diabetes")
)
prs_assocs[map, on = .(PRS), Endpoint := endpoint]

# Build set of tests
tests <- prs_assocs[,.(PRS, Target, UniProt, Gene, Endpoint)]
tests <- tests[sinfo, on = .(Target, UniProt, Gene), nomatch=0]

# Run mediation analysis using medflex package
mediation <- foreach(test_idx = tests[,.I], .combine=rbind, .export=c("soma", "hes", "prs", "pheno")) %do% {
  this_prs = tests[test_idx, PRS]
  this_target = tests[test_idx, Target]
  this_uniprot = tests[test_idx, UniProt]
  this_gene = tests[test_idx, Gene]
  this_aptvar = tests[test_idx, variable]
  this_aptamer = tests[test_idx, Aptamer]
  this_endpoint = tests[test_idx, Endpoint]

  dat <- soma[variable == this_aptvar]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]
  dat <- dat[hes[phenotype == this_endpoint], on = .(IID), nomatch=0]
  dat <- dat[pheno, on = .(IID), nomatch=0]
  dat[, prs_adj_pcs := scale(prs_adj_pcs)]
  dat[, soma_ivt_adj_batch := scale(soma_ivt_adj_batch)]
  dat[, age := scale(age)]
  dat[, sex := factor(sex)]
  dat[, event := factor(event)]

  if (dat[event == 1L, .N] < 10) return(NULL)

  extData <- neImpute(event ~ prs_adj_pcs + soma_ivt_adj_batch + age + sex, family = binomial("logit"), data=dat)
  nMod <- neModel(event ~ prs_adj_pcs0 + prs_adj_pcs1 + age + sex, family = binomial("logit"), expData=extData, se="robust")
  cf <- coef(summary(neEffdecomp(nMod)))
  ci <- confint(neEffdecomp(nMod))

  data.table(PRS = this_prs, Endpoint = this_endpoint, Target = this_target, 
						 UniProt = this_uniprot, Gene = this_gene, Aptamer = this_aptamer,
						 Mediation = rownames(cf), logOR = cf[,1], SE = cf[, 2], L95 = ci[,1],
						 U95 = ci[,2], P = cf[,4], PTE = cf[,1] / cf["total effect", 1]) 
}

# Summarise to protein level
mediation_prot <- mediation[, .(logOR = mean(logOR), SE = mean(SE), L95 = mean(L95), U95 = mean(U95), P = prot_pvalue(P, logOR), PTE = mean(PTE)),
 by = .(PRS, Endpoint, Target, UniProt, Gene, Mediation)]

# Collate
mediation <- merge(mediation_prot, mediation, by=c("PRS", "Target", "UniProt", "Gene", "Endpoint", "Mediation"), suffixes=c("", ".aptamer"))
fwrite(mediation, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/protein_mediation.txt")

# Also fit odds ratios to sanity check effect sizes
apt_or <- foreach(test_idx = tests[,.I], .combine=rbind) %dopar% {
  this_prs = tests[test_idx, PRS]
  this_target = tests[test_idx, Target]
  this_uniprot = tests[test_idx, UniProt]
  this_gene = tests[test_idx, Gene]
  this_aptvar = tests[test_idx, variable]
  this_aptamer = tests[test_idx, Aptamer]
  this_endpoint = tests[test_idx, Endpoint]

  dat <- soma[variable == this_aptvar]
  dat[hes[phenotype == this_endpoint], on = .(IID), event := event]
  dat[, soma_ivt_adj_batch := scale(soma_ivt_adj_batch)]
  dat <- dat[pheno, on = .(IID), nomatch=0]

  if (dat[,sum(event)] < 10) return(NULL)

  o1 <- glm(event ~ soma_ivt_adj_batch + age + factor(sex), data=dat, family=binomial(link="logit"))
  cf <- coef(summary(o1))
  ci <- confint(o1)

  data.table(PRS = this_prs, Endpoint = this_endpoint, Target = this_target, UniProt = this_uniprot, Gene = this_gene, Aptamer = this_aptamer,
             logOR = cf[2,1], SE = cf[2,2], OR = exp(cf[2,1]), L95 = exp(ci[2,1]), U95 = exp(ci[2,2]), P = cf[2,4])
}

# Summarise to protein level
prot_or <- apt_or[, .(logOR = mean(logOR), SE = mean(SE), OR = exp(mean(log(OR))), L95 = exp(mean(log(L95))), U95 = exp(mean(log(U95))), P = prot_pvalue(P, logOR)),
                   by = .(PRS, Endpoint, Target, UniProt, Gene)]

# Collate
prot_or <- merge(prot_or, apt_or, by=c("PRS", "Target", "UniProt", "Gene"), suffixes=c("", ".aptamer"))
fwrite(prot_or, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/protein_odds_ratios.txt")

# also get OR for PRS
tests <- unique(tests[,.(Endpoint, PRS)])
prs_or <- foreach(test_idx = tests[,.I], .combine=rbind) %dopar% {
  this_prs = tests[test_idx, PRS]
  this_endpoint = tests[test_idx, Endpoint]

  dat <- prs[PRS == this_prs & IID %in% unique(soma$IID)]
  dat[hes[phenotype == this_endpoint], on = .(IID), event := event]
  dat[, prs_adj_pcs := scale(prs_adj_pcs)]
  dat <- dat[pheno, on = .(IID), nomatch=0]

  if (dat[,sum(event)] < 10) return(NULL)

  o1 <- glm(event ~ prs_adj_pcs + age + factor(sex), data=dat, family=binomial(link="logit"))
  cf <- coef(summary(o1))
  ci <- confint(o1)

  data.table(PRS = this_prs, Endpoint = this_endpoint,
             logOR = cf[2,1], SE = cf[2,2], OR = exp(cf[2,1]), L95 = exp(ci[2,1]), U95 = exp(ci[2,2]), P = cf[2,4])
}

fwrite(prs_or, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/prs_odds_ratios.txt")

# OR adjusting for PRS
tests <- prs_assocs[,.(PRS, Target, UniProt, Gene, Endpoint)]
tests <- tests[sinfo, on = .(Target, UniProt, Gene), nomatch=0]
apt_or <- foreach(test_idx = tests[,.I], .combine=rbind) %dopar% {
  this_prs = tests[test_idx, PRS]
  this_target = tests[test_idx, Target]
  this_uniprot = tests[test_idx, UniProt]
  this_gene = tests[test_idx, Gene]
  this_aptvar = tests[test_idx, variable]
  this_aptamer = tests[test_idx, Aptamer]
  this_endpoint = tests[test_idx, Endpoint]

  dat <- soma[variable == this_aptvar]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]
  dat[hes[phenotype == this_endpoint], on = .(IID), event := event]
  dat[, soma_ivt_adj_batch := scale(soma_ivt_adj_batch)]
  dat[, prs_adj_pcs := scale(prs_adj_pcs)]

  if (dat[,sum(event)] < 10) return(NULL)

  o1 <- glm(event ~ soma_ivt_adj_batch + prs_adj_pcs, data=dat, family=binomial(link="logit"))
  cf <- coef(summary(o1))
  ci <- confint(o1)

  data.table(PRS = this_prs, Endpoint = this_endpoint, Target = this_target, UniProt = this_uniprot, Gene = this_gene, Aptamer = this_aptamer,
             coefficient = c("Protein", "PRS"), logOR = cf[-1,1], SE = cf[-1,2], OR = exp(cf[-1,1]), L95 = exp(ci[-1,1]), U95 = exp(ci[-1,2]), P = cf[-1,4])
}

# Summarise to protein level
prot_or <- apt_or[, .(logOR = mean(logOR), SE = mean(SE), OR = exp(mean(log(OR))), L95 = exp(mean(log(L95))), U95 = exp(mean(log(U95))), P = prot_pvalue(P, logOR)),
                   by = .(PRS, Endpoint, Target, UniProt, Gene, coefficient)]

# Collate
prot_or <- merge(prot_or, apt_or, by=c("PRS", "Target", "UniProt", "Gene", "coefficient"), suffixes=c("", ".aptamer"))
fwrite(prot_or, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/protein_adj_prs_odds_ratios.txt")

# Load Hazard Ratios
hes_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")
hes_assocs <- unique(hes_assocs[,.(endpoint, Gene, logHR, HR.P)])

# Build table comparing effects
med_dt <- unique(prs_assocs[,.(PRS, Endpoint, Gene, PRS.Beta)])
med_dt[prs_or, on = .(PRS, Endpoint), PRS.logOR := logOR]
med_dt <- med_dt[!is.na(PRS.logOR)]
med_dt[prot_or, on = .(PRS, Endpoint, Gene), Prot.logOR := i.logOR]
med_dt[prot_or, on = .(PRS, Endpoint, Gene), Prot.logOR.P := i.P]
med_dt[hes_assocs, on = .(Endpoint=endpoint, Gene), Prot.logHR := i.logHR]
med_dt[hes_assocs, on = .(Endpoint=endpoint, Gene), Prot.logHR.P := i.HR.P]
med_dt[mediation_prot[Mediation == "natural indirect effect"], on = .(PRS, Gene), c("Causal.effect", "Causal.P") := .(i.logOR, i.P)]
med_dt[mediation_prot[Mediation == "total effect"], on = .(PRS, Gene), PTE := Causal.effect / i.logOR]

#####################################################################################
# Test multiple mediation where multiple proteins mediate PGS to disease association
#####################################################################################

sig_med <- unique(mediation[Mediation == "natural indirect effect" & P < 0.05, .(PRS, Gene, Endpoint)])
sig_med <- sig_med[sig_med[,.N,by=.(PRS, Endpoint)][N > 1, .(PRS, Gene, Endpoint)], on = .(PRS, Gene, Endpoint)]
sig_med <- sig_med[sinfo[,.(Gene, Aptamer, variable)], on =.(Gene), nomatch=0]

multi_med <- foreach(ep = unique(sig_med$Endpoint), .combine=rbind) %do% {
  dat <- soma[variable %in% sig_med[Endpoint == ep, variable]]
  dat[, soma_ivt_adj_batch := scale(soma_ivt_adj_batch), by=variable]
  dat <- dcast(dat, IID ~ variable, value.var="soma_ivt_adj_batch")
  dat[prs[PRS == sig_med[Endpoint == ep, unique(PRS)]], on = .(IID), prs_adj_pcs := scale(i.prs_adj_pcs)]
  dat[pheno, on = .(IID), c("age", "sex") := .(scale(age), factor(sex))]
  dat[hes[phenotype == ep], on = .(IID), event := factor(i.event)]
  setDF(dat)

  m1 <- mma(x = dat[,c(sig_med[Endpoint == ep, variable], "age", "sex")],
            y = dat$event, family1=binomial(link="logit"), n2=50,
            pred = dat[,"prs_adj_pcs", drop=FALSE],
            mediator = sig_med[Endpoint == ep, variable],
            contmed = sig_med[Endpoint == ep, variable],
            jointm=list(1, sig_med[Endpoint == ep, variable]))

  indirect_effect <- m1[["a.contx"]][["estimation"]][["ie"]][1, "y1.all"]
  total_effect <- m1[["a.contx"]][["estimation"]][["te"]]["y1.prs_adj_pcs"]
  pte <- indirect_effect / total_effect

  bootstraps <- data.table(indirect_effect = m1[["a.contx"]][["bootsresults"]][["ie"]][["prs_adj_pcs"]][,"y1.all"],
	  											 total_effect = m1[["a.contx"]][["bootsresults"]][["te"]])
  bootstraps[, pte := indirect_effect / total_effect]
  bootstraps <- bootstraps[order(pte)]

  data.table(PRS = sig_med[Endpoint == ep, unique(PRS)], Endpoint = ep, PTE = pte, bootstraps)
}

fwrite(multi_med, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/mma_t2d_with_boot.txt")



