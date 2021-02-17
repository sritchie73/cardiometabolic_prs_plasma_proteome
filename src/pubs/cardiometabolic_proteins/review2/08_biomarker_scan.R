library(data.table)
library(openxlsx)
library(survival)
library(ggplot2)
library(cowplot)
library(foreach)
library(doMC)
source("src/utilities/prot_pval.R")
source("src/utilities/format_pval.R")

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

# Drop prevalent cases
soma <- soma[!prev, on = .(IID)]

# And withdrawn samples
idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
idmap <- fread(idmap_file)
soma <- soma[IID %in% idmap$Affymetrix_gwasQC_bl]

# Filter phenotype and hes data to people with somalogic proteins
pheno <- pheno[IID %in% unique(soma$IID)]
hes <- hes[IID %in% unique(soma$IID)]

# Load protein information.
sinfo <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

# Filter to aptamers passing QC
sinfo <- sinfo[Type == "Protein"]

# Select columns
sinfo <- sinfo[Type == "Protein", .(variable, SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot,
                                    Gene=Gene.Name, Entrez_id=Entrez.Gene.ID, chr, start, end, Cross_Reactivity=Characterization.Info,
                                    Mass_Spec_Confirmation=Mass.Spec.Confirmation.in.Matrix)]

# Curate aptamer sensitivity and specificity information
emilsson <- rbind(fill=TRUE,
  as.data.table(read.xlsx("data/Emilsson_etal_2018/supp_tables.xlsx", sheet="Table S3", startRow=2)),
  as.data.table(read.xlsx("data/Emilsson_etal_2018/supp_tables.xlsx", sheet="Table S4", startRow=2)))
soma_by_gene <- sinfo[,.(Gene=strsplit(Gene, "\\|")[[1]]),by=variable]
soma_by_gene[emilsson, on = .(Gene=Gene.Symbol), Mass_Spec_Confirmation := TRUE]
soma_by_gene <- soma_by_gene[!is.na(Mass_Spec_Confirmation)]

sinfo[Mass_Spec_Confirmation != "", Mass_Spec_Confirmation := "SomaLogic"]
sinfo[soma_by_gene, on = .(variable), Mass_Spec_Confirmation := ifelse(Mass_Spec_Confirmation == "", "Emilsson et al. 2018", Mass_Spec_Confirmation)]

# Are the aptamers supported by cis pQTLs?
cis_pQTLs <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
cis_pQTLs <- unique(cis_pQTLs[,.(SOMAMER_ID)])
sinfo[, cis_pQTL := ""]
sinfo[cis_pQTLs, on = .(SOMAMER_ID), cis_pQTL := "yes"]

# Fix bad entries (Aptamers for the same target with different/missing gene/uniprot information)
sinfo[Target == "14-3-3 protein family", UniProt := "P61981|Q04917"]
sinfo[Target == "Induced myeloid leukemia cell differentiation protein Mcl-1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("Q07820", "MCL1", "4170", "1", "150547027", "150552214")]
sinfo[Target == "Protein delta homolog 1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("P80370", "DLK1", "8788", "14", "101193202", "101201467")]
sinfo[Target == "Stromal cell-derived factor 1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("P48061", "CXCL12", "6387", "10", "44865601", "44880545")]

# Filter to aptamers passing QC
soma <- soma[variable %in% sinfo$variable]

###########################################
# Biomarker scan
###########################################

endpoints <- c("Myocardial infarction", "Diabetes", "Ischaemic stroke", "Atrial fibrillation")
apt_assocs <- foreach(this_endpoint = endpoints, .combine=rbind) %:% 
  foreach(this_aptvar = sinfo$variable, .combine=rbind) %dopar% {
    dat <- soma[variable == this_aptvar]
    dat <- dat[hes[phenotype == this_endpoint], on = .(IID), nomatch=0]
    dat <- dat[pheno, on = .(IID), nomatch=0]

    c1 <- coxph(Surv(followUp, event) ~ scale(soma_ivt_adj_batch) + age + factor(sex), data=dat)
    cf <- coef(summary(c1))
    ci <- confint(c1)

    data.table(endpoint = this_endpoint, Aptamer = sinfo[variable == this_aptvar, Aptamer], 
               logHR = cf[1,1], logHR.SE = cf[1,3], HR = cf[1,2], HR.L95 = exp(ci[1,1]), 
               HR.U95 = exp(ci[1,2]), HR.P = cf[1,5])
}

# Add in information about the aptamer target
apt_assocs <- sinfo[apt_assocs, on = .(Aptamer)]

# Average associations across aptamers for each protein
prot_assocs <- apt_assocs[, .(logHR = mean(logHR), logHR.SE = mean(logHR.SE), 
  HR = exp(mean(log(HR))), HR.L95 = exp(mean(log(HR.L95))), HR.U95 = exp(mean(log(HR.U95))),
  HR.P = prot_pvalue(HR.P, logHR)),
  by=.(endpoint, Target, UniProt, Gene)]
prot_assocs[, HR.FDR := p.adjust(HR.P, method="fdr"), by=.(endpoint)]

# Merge protein and aptamer associations to a single table
prot_assocs <- merge(prot_assocs, apt_assocs, by = c("endpoint", "Target", "UniProt", "Gene"), suffixes=c("", ".aptamer"))

# Select columns
prot_assocs <- prot_assocs[, .(endpoint, Target, UniProt, Gene, Entrez_id, chr, start, end,
  logHR, logHR.SE, HR, HR.L95, HR.U95, HR.P, HR.FDR, Aptamer, logHR.aptamer, logHR.SE.aptamer, 
  HR.aptamer, HR.L95.aptamer, HR.U95.aptamer, HR.P.aptamer, Cross_Reactivity, Mass_Spec_Confirmation,
  cis_pQTL)]

# Order rows
prot_assocs <- prot_assocs[order(abs(logHR.aptamer))][order(HR.P.aptamer)][order(abs(logHR))][order(HR.P)][order(HR.FDR)][order(endpoint)]

# Write out
fwrite(prot_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")

###########################################
# Compare overlap with PRS associations
###########################################

prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- prs_assocs[, .(PRS, Target, UniProt, Gene, Beta, L95, U95, P, FDR)]
prs_assocs <- unique(prs_assocs)

prot_bio <- prot_assocs[,.(endpoint, Target, UniProt, Gene, HR, HR.L95, HR.U95, HR.P, HR.FDR)]
prot_bio <- unique(sig_bio)

prs_endmap <- data.table(
  endpoint = c("Atrial fibrillation", "Myocardial infarction", "Ischaemic stroke", "Diabetes"),
  PRS = c("AF_PRS", "CAD_PRS", "IS_PRS", "T2D_PRS")
)

prs_assocs[prs_endmap, on = .(PRS), endpoint := i.endpoint]
prot_bio[prs_endmap, on = .(endpoint), PRS := i.PRS]

comp <- merge(prs_assocs, prot_bio, by=c("PRS", "endpoint", "Target", "UniProt", "Gene"))
comp <- comp[FDR < 0.05 | HR.FDR < 0.05]
comp[FDR < 0.05 & HR.FDR < 0.05, anno := "both"]
comp[FDR < 0.05 & HR.FDR >= 0.05, anno := "PRS"]
comp[FDR >= 0.05 & HR.FDR < 0.05, anno := "endpoint"]

comp[, discordant := sign(log(HR)) != sign(Beta)]

# Get correlation labels
tests <- unique(comp[,.(PRS, endpoint)])
cor_assocs <- foreach(test_idx = tests[,.I], .combine=rbind) %do% {
  this_prs <- tests[test_idx, PRS]
  this_end <- tests[test_idx, endpoint]

  dat <- comp[PRS == this_prs & endpoint == this_end]
  if (nrow(dat) > 1) {
    c1 <- dat[,cor.test(Beta, HR)]
    data.table(PRS = this_prs, endpoint = this_end, r = c1$estimate, 
               L95 = c1$conf.int[1], U95 = c1$conf.int[2], P = c1$p.value)
  }
}
cor_assocs[, label := sprintf("Pearson r = %s\n[95%% CI: %s-%s]\nP=%s",
             round(r*100)/100, round(L95*100)/100, round(U95*100)/100, 
             my_format_pval(P))]

g <- ggplot(comp) +
  aes(x = Beta, xmin = L95, xmax = U95,
      y = HR, ymin = HR.L95, ymax = HR.U95, 
      color = anno) +
	geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
	geom_errorbarh(height=0, alpha=0.8, size=0.5) +
	geom_errorbar(width=0, alpha=0.8, size=0.5) +
	geom_point(shape = 19, size=1.3) +
	geom_text(data=cor_assocs, inherit.aes=FALSE, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1, aes(label=label)) +
  facet_wrap(~ PRS + endpoint) +
	scale_x_continuous(name = "PRS Beta (95% CI)") +
	scale_y_log10(name = "Protein HR (95% CI)") +
	theme_bw() + theme(
		axis.title=element_text(size=8), axis.text=element_text(size=8),
		panel.grid=element_blank(), legend.position="bottom"
	)
ggsave(g, width=5.4, height=2.8, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/PRS_HR_compare.pdf")
  
  




                       
                               
