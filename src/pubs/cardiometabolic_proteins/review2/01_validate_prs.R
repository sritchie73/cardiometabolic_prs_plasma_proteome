library(data.table)
library(survival)
library(ggplot2)
library(cowplot)

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
  CAD_metaGRS = fread("analyses/GRS_profiles/CAD_metaGRS/profile.sscore.gz"),
  Stroke_metaGRS = fread("analyses/GRS_profiles/Stroke_metaGRS/profile.sscore.gz"),
  T2D_PRS = fread("analyses/GRS_profiles/T2D_2018/profile.sscore.gz"),
  Afib_PRS = fread("analyses/GRS_profiles/Afib_2018/profile.sscore.gz"),
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

# Load phenotype data - use the HES data because it has all samples mapped
pheno_file <- list.files("data/INTERVAL/HES", pattern="INTERVALdata_[^(p|P)3].*", full.names=TRUE)
pheno <- fread(pheno_file)

idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
idmap <- fread(idmap_file)

pheno[idmap, on = .(identifier), IID := Affymetrix_gwasQC_bl]
pheno <- pheno[!is.na(IID), .(IID, age=agePulse, sex=sexPulse)]

pheno <- pheno[!prev, on = .(IID)]

# Build wide table for model fitting:
dat <- dcast(prs, IID ~ PRS, value.var="prs_adj_pcs")
dat <- dat[pheno, on = .(IID), nomatch=0]
dat <- dat[IID %in% unique(hes$IID)]

dat[hes[phenotype == "Myocardial infarction"], on = .(IID), c("MI_event", "MI_follow") := .(event, followUp)]
dat[hes[phenotype == "Diabetes"], on = .(IID), c("Diab_event", "Diab_follow") := .(event, followUp)]
dat[hes[phenotype == "Ischaemic stroke"], on = .(IID), c("IS_event", "IS_follow") := .(event, followUp)]
dat[hes[phenotype == "Atrial fibrillation"], on = .(IID), c("Afib_event", "Afib_follow") := .(event, followUp)]
dat[hes[phenotype == "End stage renal disease"], on = .(IID), c("ESRD_event", "ESRD_follow") := .(event, followUp)]

# Load creatinine and eGFR information
source("src/pubs/cardiometabolic_proteins/review2/compute_eGFR.R")
dat <- merge(dat, crea, by = "IID", all.x=TRUE)

# Flag which samples have somalogic proteins 
soma <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")
dat[, has_soma := FALSE]
dat[soma, on = .(IID), has_soma := TRUE]

# Fit cox models for endpoint and extract coefficients
cox.test <- function(formula, data) {
  cx <- coxph(formula, data)
  cf <- coef(summary(cx))
  ci <- confint(cx)
  data.table(logHR=cf[1,1], logHR.SE=cf[1,3], HR=cf[1,2], HR.L95=exp(ci[1,1]), HR.U95=exp(ci[1,2]), HR.P=cf[1,5])
}

prs_hrs <- rbind(idcol="Cohort",
  "Full"=rbind(idcol="test",
    "MI ~ CAD_PRS" = cbind(cox.test(Surv(MI_follow, MI_event) ~ scale(CAD_metaGRS) + factor(sex) + age, data=dat), samples=dat[,.N], events=dat[,sum(MI_event)]),
    "IS ~ Stroke_PRS" = cbind(cox.test(Surv(IS_follow, IS_event) ~ scale(Stroke_metaGRS) + factor(sex) + age, data=dat), samples=dat[,.N], events=dat[,sum(IS_event)]),
    "Diab ~ T2D_PRS" = cbind(cox.test(Surv(Diab_follow, Diab_event) ~ scale(T2D_PRS) + factor(sex) + age, data=dat), samples=dat[,.N], events=dat[,sum(Diab_event)]),
    "Afib ~ Afib_PRS" = cbind(cox.test(Surv(Afib_follow, Afib_event) ~ scale(Afib_PRS) + factor(sex) + age, data=dat), samples=dat[,.N], events=dat[,sum(Afib_event)]),
    "ESRD ~ CKD_PRS" = cbind(cox.test(Surv(ESRD_follow, ESRD_event) ~ scale(CKD_PRS) + factor(sex) + age, data=dat), samples=dat[,.N], events=dat[,sum(ESRD_event)])
  ),
  "Soma"=rbind(idcol="test",
    "MI ~ CAD_PRS" = cbind(cox.test(Surv(MI_follow, MI_event) ~ scale(CAD_metaGRS) + factor(sex) + age, data=dat[(has_soma)]), samples=dat[(has_soma),.N], events=dat[(has_soma),sum(MI_event)]),
    "IS ~ Stroke_PRS" = cbind(cox.test(Surv(IS_follow, IS_event) ~ scale(Stroke_metaGRS) + factor(sex) + age, data=dat[(has_soma)]), samples=dat[(has_soma),.N], events=dat[(has_soma),sum(IS_event)]),
    "Diab ~ T2D_PRS" = cbind(cox.test(Surv(Diab_follow, Diab_event) ~ scale(T2D_PRS) + factor(sex) + age, data=dat[(has_soma)]), samples=dat[(has_soma),.N], events=dat[(has_soma),sum(Diab_event)]),
    "Afib ~ Afib_PRS" = cbind(cox.test(Surv(Afib_follow, Afib_event) ~ scale(Afib_PRS) + factor(sex) + age, data=dat[(has_soma)]), samples=dat[(has_soma),.N], events=dat[(has_soma),sum(Afib_event)]),
    "ESRD ~ CKD_PRS" = cbind(cox.test(Surv(ESRD_follow, ESRD_event) ~ scale(CKD_PRS) + factor(sex) + age, data=dat[(has_soma)]), samples=dat[(has_soma),.N], events=dat[(has_soma),sum(ESRD_event)])
  )
)

# Test the CKD PRS against eGFR
lm.test <- function(formula, data) {
  l1 <- lm(formula, data)
  cf <- coef(summary(l1))
  ci <- confint(l1)
  data.table(beta=cf[2,1], se=cf[2,2], L95=ci[2,1], U95=ci[2,2], P=cf[2,4])
}

ckd_egfr <- rbind(idcol="Cohort",
  "Full"=rbind(idcol="test",
     "eGFR ~ CKD PRS"=cbind(lm.test(eGFR ~ scale(CKD_PRS) + age + factor(sex), data=dat[!is.na(eGFR)]), samples=dat[!is.na(eGFR),.N])
   ),
  "Soma"=rbind(idcol="test",
     "eGFR ~ CKD PRS"=cbind(lm.test(eGFR ~ scale(CKD_PRS) + age + factor(sex), data=dat[!is.na(eGFR) & (has_soma)]), samples=dat[!is.na(eGFR) & (has_soma),.N])
   )
)

# Write out information
system("mkdir -p analyses/pub/cardiometabolic_proteins/review2", wait=TRUE)
fwrite(prs_hrs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/prs_hes_associations.txt")
fwrite(ckd_egfr, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/ckd_prs_egfr_associations.txt")

# Create plots
g1 <- ggplot(prs_hrs) + 
  aes(x=HR, xmin=HR.L95, xmax=HR.U95, y=factor(test)) +
  geom_vline(xintercept=1, linetype=2) + 
  geom_errorbarh(height=0, alpha=0.8) +
  geom_point(shape=18, size=3) + 
  facet_wrap(~ factor(Cohort, levels=c("Soma", "Full")), ncol=2) +
  xlab("Hazard Ratio per SD increase of PRS") +
  ylab("") + 
  theme_bw()

g2 <- ggplot(ckd_egfr) + 
  aes(x=beta, xmin=L95, xmax=U95, y=factor(test)) +
  geom_vline(xintercept=0, linetype=2) + 
  geom_errorbarh(height=0, alpha=0.8) +
  geom_point(shape=18, size=3) + 
  facet_wrap(~ factor(Cohort, levels=c("Soma", "Full")), ncol=2) +
  xlab("CKD-EPI eGFR (mL/min/1.73 m2) per SD of PRS") +
  ylab("") + 
  theme_bw()

g <- plot_grid(g1, g2, nrow=2, rel_heights=c(0.7, 0.3))
ggsave(g, width=6.2, height=5, file="analyses/pub/cardiometabolic_proteins/review2/prs_validation.pdf")




