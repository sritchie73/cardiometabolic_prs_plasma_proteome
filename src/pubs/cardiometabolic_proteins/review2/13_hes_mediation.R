library(data.table)
library(foreach)
library(doMC)
library(survival)
library(mediation)
library(powerMediation)
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
pheno <- pheno[!is.na(IID), .(IID, age=agePulse, sex=sexPulse)]

pheno <- pheno[!prev, on = .(IID)]

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

# Load significant PRS to protein associations
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[FDR < 0.05,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])

# Filter to aptamers of interest
sinfo <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
sinfo <- sinfo[, .(variable, SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, Gene=Gene.Name)]
sinfo <- sinfo[Target %in% prs_assocs$Target]
soma <- soma[variable %in% sinfo$variable]

# Load in PRS to endpoint associations
prs_hes <- fread("analyses/pub/cardiometabolic_proteins/review2/prs_hes_associations.txt")
prs_hes <- prs_hes[Cohort == "Soma"]

# Map between tables
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  endpoint = c("Atrial fibrillation", "Myocardial infarction", "End stage renal disease", "Ischaemic stroke", "Diabetes"),
  prs_hes_test = c("Afib ~ Afib_PRS", "MI ~ CAD_PRS", "ESRD ~ CKD_PRS", "IS ~ Stroke_PRS", "Diab ~ T2D_PRS")
)
prs_hes[map, on = .(test=prs_hes_test), c("PRS", "Endpoint") := .(PRS, endpoint)]

# Drop CKD (no endpoints) and AF (no associated proteins) from PRS to endpoint associations table
prs_hes <- prs_hes[PRS %in% prs_assocs$PRS]
prs_hes <- prs_hes[!is.na(logHR)]

# Build set of tests
tests <- prs_assocs[,.(PRS, Target, UniProt, Gene)]
tests <- tests[sinfo, on = .(Target, UniProt, Gene)]
tests[map, on = .(PRS), Endpoint := i.endpoint]
tests <- tests[prs_hes[,.(PRS, Endpoint)], on = .(PRS, Endpoint), nomatch=0]

############################ --- Finished loading data

# Compute proportion of variance in PRS explained by protein
# using mediation analysis package (equivalent to proportion of
# treatment effect, e.g. 1 - logHR.adj / logHR.
mediation <- foreach(test_idx = tests[,.I], .combine=rbind) %dopar% {
  this_prs = tests[test_idx, PRS]
  this_target = tests[test_idx, Target]
  this_uniprot = tests[test_idx, UniProt]
  this_gene = tests[test_idx, Gene]
  this_aptvar = tests[test_idx, variable]
  this_aptamer = tests[test_idx, Aptamer]
  this_endpoint = tests[test_idx, Endpoint]

  dat <- soma[variable == this_aptvar]
  dat <- dat[pheno, on = .(IID), nomatch=0]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]
  dat[hes[phenotype == this_endpoint], on = .(IID), c("event", "followUp") := .(event, followUp)]
  dat[, prs_adj_pcs := scale(prs_adj_pcs)]
  dat[, soma_ivt_adj_batch := scale(soma_ivt_adj_batch)]

  setDF(dat)
  l1 <- lm(soma_ivt_adj_batch ~ prs_adj_pcs, data=dat) 
  cx <- survreg(Surv(followUp, event) ~ prs_adj_pcs + soma_ivt_adj_batch + age + factor(sex), data=dat) 

  m1 <- mediate(l1, cx, treat="prs_adj_pcs", mediator="soma_ivt_adj_batch")
  cf <- summary(m1)

  data.table(PRS = this_prs, Target = this_target, UniProt = this_uniprot, Gene = this_gene, Aptamer = this_aptamer,
             PRS.effect = cf$z.avg, PRS.effect.L95 = cf$z.avg.ci[1], PRS.effect.U95 = cf$z.avg.ci[2], PRS.effect.P = cf$z.avg.p,
             Prot.effect = cf$z.avg, Prot.effect.L95 = cf$z.avg.ci[1], Prot.effect.U95 = cf$z.avg.ci[2], Prot.effect.P = cf$z.avg.p,
             PTE = cf$n.avg, PTE.L95 = cf$n.avg.ci[1], PTE.U95 = cf$n.avg.ci[2], PTE.P = cf$n.avg.p)
}

# Summarise to protein level
prs_prot_pte <- mediation[, .(PTE = mean(PTE), PTE.L95 = mean(PTE.L95), PTE.U95 = mean(PTE.U95),
                              PTE.P = prot_pvalue(PTE.P, PTE)),
                          by = .(PRS, Target, UniProt, Gene)]

# Plot comparison of protein to disease associations before and after adjustment
comp <- prot_med[coefficient != PRS,.(Endpoint, PRS, Target, UniProt, Gene, HR, HR.L95, HR.U95, HR.P)]
comp <- merge(hes_assocs, comp, by=c("Endpoint", "PRS", "Target", "UniProt", "Gene"), suffixes=c("", ".adj"))

# plot limits
rdt <- comp[, .(x=c(min(c(HR.L95, HR.L95.adj)), 1,
									  max(c(HR.U95, HR.U95.adj))),
									  ymax=c(min(c(HR.L95, HR.L95.adj)), 1,
										 			 max(c(HR.U95, HR.U95.adj))),
									  ymin=c(1,1,1))]

# Expand
rdt[, c("x", "ymax") := .(log(x), log(ymax))]
rdt[, x := x + sign(x)*sum(abs(x))*0.05]
rdt[, ymax := ymax + sign(ymax)*sum(abs(ymax))*0.05]
rdt[, c("x", "ymax") := .(exp(x), exp(ymax))]

# Plot
g <- ggplot(comp) +
  aes(x=HR, xmin=HR.L95, xmax=HR.U95, 
      y=HR.adj, ymin=HR.L95.adj, ymax=HR.U95.adj) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
  geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.7, size=0.5) +
  geom_errorbar(width=0, alpha=0.7, size=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#737373") +
  geom_point(shape = 19) +
  facet_wrap(~ PRS) +
  scale_x_log10(name = "HR (95% CI)", expand=c(0,0)) +
  scale_y_log10(name = "HR (95% CI)", expand=c(0,0)) +
  theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank()
  )
ggsave(g, width=5.1, height=2.1, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/prot_hes_prs_adjusted.pdf")

# Build wide table for supp
supp <- merge(prot_med[coefficient != PRS], prot_med[coefficient == PRS], 
              by = c("Endpoint", "PRS", "Target", "UniProt", "Gene"), 
              suffixes=c(".protein", ".PRS"))

# Compute proportion of variance in PRS explained by all
# PRS-associated proteins
all_med <- foreach(this_prs = unique(tests$PRS), .combine=rbind) %dopar% {
  this_targets = tests[PRS == this_prs, Target]
  this_uniprots = tests[PRS == this_prs, UniProt]
  this_genes = tests[PRS == this_prs, Gene]
  this_aptvars = tests[PRS == this_prs, variable]
  this_aptamers = tests[PRS == this_prs, Aptamer]
  this_endpoint = tests[PRS == this_prs, unique(Endpoint)]
  if (length(unique(this_genes)) == 1) {
    return(NULL)
  }

  dat <- soma[variable %in% this_aptvars]
  dat[, soma_ivt_adj_batch := scale(soma_ivt_adj_batch), by=variable]
  dat <- dcast(dat, IID ~ variable, value.var="soma_ivt_adj_batch")
  dat <- dat[pheno, on = .(IID), nomatch=0]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]
  dat[, prs_adj_pcs := scale(prs_adj_pcs)]
  dat[hes[phenotype == this_endpoint], on = .(IID), c("event", "followUp") := .(event, followUp)]
  dat[, sex := sex - 1]
  setDF(dat)

  m1 <- mma(x = dat[,c(this_aptvars, "age", "sex")], 
            y = Surv(dat$followUp, dat$event), 
            pred = dat[,"prs_adj_pcs", drop=FALSE], 
            type = "lp", n = 40, n2 = 40,
            mediator = this_aptvars, 
            contmed = this_aptvars, 
            jointm=list(1, this_aptvars))

  indirect_effect <- m1[["a.contx"]][["estimation"]][["ie"]][1, "y1.all"]
  total_effect <- m1[["a.contx"]][["estimation"]][["te"]]["y1.prs_adj_pcs"]
  pte <- indirect_effect / total_effect

  bootstraps <- data.table(indirect_effect = m1[["a.contx"]][["bootsresults"]][["ie"]][["prs_adj_pcs"]][,"y1.all"],
                           total_effect = m1[["a.contx"]][["bootsresults"]][["te"]])
  bootstraps[, pte := indirect_effect / total_effect]
  bootstraps <- bootstraps[order(pte)]
  ci95 <- bootstraps[-c(1, .N)][c(1, .N), pte]

  data.table(PRS = this_prs, PTE = pte, PTE.L95 = ci95[1], PTE.U95 = ci95[2])
}

# Get P-value from 95% CI
# https://www.bmj.com/content/343/bmj.d2304
p_from_ci <- function(est, l95, u95) {
  se = (u95 - l95) / (2 * 1.96)
  z = abs(est/se)
  p = exp(-0.717 * z - 0.416 * z^2)
  return(p)
}
all_med[, PTE.P := p_from_ci(PTE, PTE.L95, PTE.U95)]

prs_prot_pte <- rbind(idcol="type", fill=TRUE, single=prs_prot_pte, all=all_med)
prs_prot_pte[type == "all", c("Target", "UniProt", "Gene") := "All"]
prs_prot_pte <- prs_prot_pte[order(PTE)][order(-type)][order(PRS)]
prs_prot_pte[,plot_order := .I]

# Plot
g <- ggplot(prs_prot_pte) +
  aes(x=factor(plot_order), y=PTE, ymin=PTE.L95, ymax=PTE.U95, color=type) +
  geom_hline(yintercept=0, linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_y_continuous(name="PTE (95% CI)", limits=c(-0.10, 0.50), oob=oob_keep) +
  scale_x_discrete(name="", labels=prs_prot_pte$Gene, breaks=prs_prot_pte$plot_order) +
  facet_grid(. ~ PRS, space="free_x", scales="free", shrink=TRUE) +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        legend.position="bottom")
ggsave(g, width=7.2, height=3, file="analyses/pub/cardiometabolic_proteins/review2/mediation.pdf", useDingbats=FALSE)

# Quantify PRS to protein associations after adjusting for all proteins in Cox regression
all_med_cox <- foreach(this_prs = unique(tests$PRS), .combine=rbind) %dopar% {
  this_targets = tests[PRS == this_prs, Target]
  this_uniprots = tests[PRS == this_prs, UniProt]
  this_genes = tests[PRS == this_prs, Gene]
  this_aptvars = tests[PRS == this_prs, variable]
  this_aptamers = tests[PRS == this_prs, Aptamer]
  this_endpoint = tests[PRS == this_prs, unique(Endpoint)]
  if (length(unique(this_genes)) == 1) {
    return(NULL)
  }

  dat <- soma[variable %in% this_aptvars]
  dat <- dcast(dat, IID ~ variable, value.var="soma_ivt_adj_batch")
  dat <- dat[pheno, on = .(IID), nomatch=0]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]
  dat[hes[phenotype == this_endpoint], on = .(IID), c("event", "followUp") := .(event, followUp)]

  mf <- sprintf("Surv(followUp, event) ~ scale(prs_adj_pcs) + %s + age + factor(sex)",
                paste(sprintf("scale(%s)", this_aptvars), collapse=" + "))

  cx <- coxph(as.formula(mf), data=dat)
  cf <- coef(summary(cx))
  ci <- confint(cx)

  nvars <- length(this_aptvars) + 1
  cf <- data.table(coefficient = c(this_prs, this_aptvars), logHR = cf[1:nvars, 1], logHR.SE = cf[1:nvars, 3],
                   HR = cf[1:nvars, 2], HR.L95 = exp(ci[1:nvars, 1]), HR.U95 = exp(ci[1:nvars, 2]), HR.P = cf[1:nvars, 5])
  info <- data.table(Endpoint = this_endpoint, PRS = this_prs, Target = c("All", this_targets), 
                     UniProt = c("All", this_uniprots), Gene = c("All", this_genes))
  cbind(info, cf)
}

prs_adj_all <- all_med_cox[coefficient == PRS, .(Endpoint, PRS, Target, UniProt, Gene,
                           logHR.PRS=logHR, logHR.SE.PRS=logHR.SE, HR.PRS=HR, HR.L95.PRS=HR.L95, 
                           HR.U95.PRS=HR.U95, HR.P.PRS=HR.P)]
supp <- rbind(supp, prs_adj_all, fill=TRUE)

# Add PTE to supp table
supp[prs_hes, on = .(PRS), mediated_effect := i.logHR]
supp[prs_prot_pte, on = .(PRS, UniProt, Target, Gene), c("PTE", "PTE.L95", "PTE.U95", "PTE.P") := .(PTE, PTE.L95, PTE.U95, PTE.P)]
supp[, mediated_effect := mediated_effect * PTE]

# Compute power for each test
supp[prs_assocs, on = .(PRS, Target), PRS.Protein.Beta := i.PRS.Beta]
supp[, ncases := fcase(PRS == "CAD_PRS", 15, PRS == "IS_PRS", 3, PRS == "T2D_PRS", 27)]

# For multiple- mediation need to know correlation between PRS and all aptamers
prs_all_r <- foreach(this_prs = unique(tests$PRS), .combine=rbind) %do% {
  this_targets = tests[PRS == this_prs, Target]
  this_uniprots = tests[PRS == this_prs, UniProt]
  this_genes = tests[PRS == this_prs, Gene]
  this_aptvars = tests[PRS == this_prs, variable]
  this_aptamers = tests[PRS == this_prs, Aptamer]
  this_endpoint = tests[PRS == this_prs, unique(Endpoint)]
  if (length(unique(this_genes)) == 1) {
    return(NULL)
  }

  dat <- soma[variable %in% this_aptvars]
  dat <- dcast(dat, IID ~ variable, value.var="soma_ivt_adj_batch")
  dat <- dat[pheno, on = .(IID), nomatch=0]
  dat <- dat[prs[PRS == this_prs], on = .(IID), nomatch=0]

  mf <- sprintf("scale(prs_adj_pcs) ~ %s", paste(sprintf("scale(%s)", this_aptvars), collapse=" + "))
  l1 <- lm(as.formula(mf), data=dat)

  data.table(PRS = this_prs, Target="All", r2 = summary(l1)$r.squared, r = sqrt(summary(l1)$r.squared))
}
supp[prs_all_r, on = .(PRS, Target), PRS.Protein.Beta := i.r]

supp[, PTE.power := powerMediation.VSMc.cox(
       n=3087, b2 = mediated_effect, sigma.m = 1, psi = ncases/3087,
       corr.xm = PRS.Protein.Beta, verbose = FALSE
     )$power]
supp[, c("PRS.Protein.Beta", "ncases") := NULL]

# Write out
fwrite(supp, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/prs_hes_prot_mediation.txt")

# Compute power for a range of scenarios and sample sizes
sample_sizes <- expand.grid(1:10, 3:5)
sample_sizes <- sample_sizes$Var1 * 10^(sample_sizes$Var2)
prevalence <- c(0.01, 0.15, 0.5)
delta_loghr <- c(.01, .04, .1)

power <- foreach(ii = sample_sizes, .combine=rbind) %:%
  foreach(jj = prevalence, .combine=rbind) %:%
    foreach(kk = delta_loghr, .combine=rbind) %do% {
      res <- data.table(n_samples = ii, prevalence = jj, delta_loghr = kk)
      res[, power := powerMediation.VSMc.cox(
        n = n_samples, b2 = delta_loghr, sigma.m = 1, psi = prevalence,
        corr.xm =  0.076, verbose = FALSE
      )$power]
}

power_multi_med <- foreach(ii = sample_sizes, .combine=rbind) %:%
  foreach(jj = prevalence, .combine=rbind) %do% {
      res <- data.table(n_samples = ii, prevalence = jj, delta_loghr = 0.295)
      res[, power := powerMediation.VSMc.cox(
        n = n_samples, b2 = delta_loghr, sigma.m = 1, psi = prevalence,
        corr.xm =  0.2, verbose = FALSE
      )$power]
}
power <- rbind(power, power_multi_med)

g <- ggplot(power) +
  aes(x=n_samples, y=power, color = factor(delta_loghr), linetype = factor(prevalence)) +
  geom_vline(xintercept=sample_sizes, color="#ebebeb", size=0.3) +
  geom_hline(yintercept=0.8, color="#252525", linetype=2) +
  geom_line() +
  scale_color_manual(values=c("0.01"="#e41a1c", "0.04"="#377eb8", "0.1"="#4daf4a", "0.295"="#ff6600")) +
  scale_x_log10(name="Sample size", breaks=10^(3:6), labels=c("1K", "10K", "100K", "1M")) +
  scale_y_continuous(name="Power", breaks=seq(0, 1, by=0.2)) +
  theme_bw() +
  theme(legend.position="bottom", panel.grid.minor=element_blank(),
        axis.title=element_text(size=10), axis.text=element_text(size=8),
        legend.text=element_text(size=8), legend.title=element_blank())
ggsave(g, width=2.1, height=2.2, file="analyses/pub/cardiometabolic_proteins/review2/prs_mediation_power.pdf")



