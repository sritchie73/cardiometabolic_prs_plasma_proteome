library(data.table)
library(foreach)
library(ggplot2)
source("src/utilities/format_pval.R")

# Load instruments and coloc results
ivs <- fread("analyses/pub/cardiometabolic_proteins/review2/cis_ivs.txt")
coloc <- fread("analyses/pub/cardiometabolic_proteins/review2/cis_iv_coloc.txt")

ivs[coloc, on = .(PRS, GWAS, Aptamer, pQTL.chr, pQTL.pos), 
  c("coloc.nsnps", "coloc.PP.H0", "coloc.PP.H1", "coloc.PP.H2", "coloc.PP.H3", "coloc.PP.H4", "colocalizes") :=
  .(nsnps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, colocalises)]

# Load and collate MR results
mr_apt <- fread("analyses/pub/cardiometabolic_proteins/review2/mr_all_aptamers.txt")
mr_prot <- fread("analyses/pub/cardiometabolic_proteins/review2/mr_all_proteins.txt")
mr_dt <- merge(mr_prot, mr_apt, by = c("PRS", "GWAS", "Target", "UniProt", "Gene", "Method"), suffixes=c("", ".aptamer")) 

# Load MR summary
mr_summary <- fread("analyses/pub/cardiometabolic_proteins/review2/mr_summary_all_proteins.txt")
mr_summary[ivs[(colocalizes)], on = .(PRS, GWAS, Target, UniProt, Gene), colocalization := TRUE]
mr_summary[is.na(colocalization), colocalization := FALSE]
mr_dt <- merge(mr_summary, mr_dt, by = c("PRS", "GWAS", "Target", "UniProt", "Gene"), suffixes=c("", ".protein"))

# Order rows 
mr_dt[, Method := factor(Method, levels=c("IVW", "Median", "Weighted Median", "Weighted Mode", "MR Egger", "(Intercept)"))]
mr_dt <- mr_dt[order(Method)][order(mr_pval)][order(mr_fdr)][order(pleiotropy_pval < 0.05)][order(PRS)]
ivs <- ivs[order(pQTL.P)]
ivs <- ivs[unique(mr_dt[,.(PRS, GWAS, Target, Gene, UniProt, Aptamer)]), on = .(PRS, GWAS, Target, Gene, UniProt, Aptamer)]

# Write out collated tables for supp
fwrite(ivs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/collated_ivs.txt")
fwrite(mr_dt, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/collated_mr.txt")
 
#####################################
# Compare to PRS and HES assocs
#####################################

# Load prs and hes assocations
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
hes_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")

# Filter to protein associations
prs_assocs <- unique(prs_assocs[,.(PRS, Target, UniProt, Gene, Beta, L95, U95, P, FDR)])
hes_assocs <- unique(hes_assocs[,.(endpoint, Target, UniProt, Gene, HR, HR.L95, HR.U95, HR.P, HR.FDR)])
mr_assocs <- unique(mr_dt[,.(GWAS, Target, UniProt, Gene, MR.OR=exp(mr_estimate), MR.OR.L95=exp(mr_L95), MR.OR.U95=exp(mr_U95), 
                             MR.P=mr_pval, MR.FDR=mr_fdr, MR.pleiotropy=pleiotropy_pval, colocalization)])

# Map between outcomes
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  endpoint = c("Atrial fibrillation", "Myocardial infarction", "End stage renal disease", "Ischaemic stroke", "Diabetes"),
  GWAS = c("Afib", "CAD", "CKD", "StrokeIS", "T2DadjBMI")
)

prs_assocs <- map[prs_assocs, on = .(PRS)]
hes_assocs <- map[hes_assocs, on = .(endpoint), nomatch=0]
mr_assocs <- map[mr_assocs, on = .(GWAS)]

# Get set of FDR significant associations across all tests
sig_set <- unique(rbind(
  prs_assocs[FDR < 0.05, .(PRS, endpoint, GWAS, Target, UniProt, Gene)],
  hes_assocs[HR.FDR < 0.05, .(PRS, endpoint, GWAS, Target, UniProt, Gene)],
  mr_assocs[MR.FDR < 0.05 & MR.pleiotropy > 0.05, .(PRS, endpoint, GWAS, Target, UniProt, Gene)]
))

comp <- mr_assocs[sig_set, on = .(PRS, endpoint, GWAS, Target, UniProt, Gene), nomatch=0]
comp <- merge(comp, prs_assocs, by=c("PRS", "endpoint", "GWAS", "Target", "UniProt", "Gene"), all.x=TRUE)
comp <- merge(comp, hes_assocs, by=c("PRS", "endpoint", "GWAS", "Target", "UniProt", "Gene"), all.x=TRUE)

# Flag discordant directions
comp[, discordant.hes.prs := sign(Beta) != sign(log(HR))]
comp[, discordant.mr.prs := sign(Beta) != sign(log(MR.OR))]
comp[, discordant.mr.hes := sign(log(MR.OR)) != sign(log(HR))]
comp[, concordant := !(discordant.hes.prs) & !(discordant.mr.prs) & !(discordant.mr.hes)]

# Plot MR vs PRS
comp[, anno := "drop"]
comp[FDR < 0.05 & MR.FDR < 0.05 & MR.pleiotropy > 0.05, anno := "both"]
comp[FDR < 0.05 & MR.FDR >= 0.05, anno := "PRS"]
comp[FDR >= 0.05 & MR.FDR < 0.05 & MR.pleiotropy > 0.05, anno := "MR"]

tests <- unique(comp[anno != "drop",.(PRS, GWAS)])
cor_assocs <- foreach(test_idx = tests[,.I], .combine=rbind) %do% {
  this_prs <- tests[test_idx, PRS]
  this_gwas <- tests[test_idx, GWAS]

  dat <- comp[PRS == this_prs & GWAS == this_gwas]
  if (nrow(dat) > 2) {
    c1 <- dat[,cor.test(Beta, MR.OR)]
    data.table(PRS = this_prs, GWAS = this_gwas, r = c1$estimate,
               L95 = c1$conf.int[1], U95 = c1$conf.int[2], P = c1$p.value)
  }
}
cor_assocs[, label := sprintf("Pearson r = %s\n[95%% CI: %s-%s]\nP=%s",
             round(r*100)/100, round(L95*100)/100, round(U95*100)/100,
             my_format_pval(P))]

g <- ggplot(comp[anno != "drop"]) +
  aes(x = Beta, xmin = L95, xmax = U95,
      y = MR.OR, ymin = MR.OR.L95, ymax = MR.OR.U95,
      color = anno) +
  geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.8, size=0.5) +
  geom_errorbar(width=0, alpha=0.8, size=0.5) +
  geom_point(shape = 19, size=1.3) +
  geom_text(data=cor_assocs, inherit.aes=FALSE, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1, aes(label=label)) +
  facet_wrap(~ PRS, nrow=1) +
  scale_x_continuous(name = "PRS Beta (95% CI)") +
  scale_y_log10(name = "Causal OR (95% CI)") +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g, width=7.2, height=2.4, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/PRS_MR_compare.pdf")

# Plot MR vs HES
comp[, anno := "drop"]
comp[HR.FDR < 0.05 & MR.FDR < 0.05 & MR.pleiotropy > 0.05, anno := "both"]
comp[HR.FDR < 0.05 & MR.FDR >= 0.05, anno := "HES"]
comp[HR.FDR >= 0.05 & MR.FDR < 0.05 & MR.pleiotropy > 0.05, anno := "MR"]

tests <- unique(comp[anno != "drop",.(PRS, endpoint)])
cor_assocs <- foreach(test_idx = tests[,.I], .combine=rbind) %do% {
  this_prs <- tests[test_idx, PRS]
  this_endpoint <- tests[test_idx, endpoint]

  dat <- comp[PRS == this_prs & endpoint == this_endpoint]
  if (nrow(dat) > 2) {
    c1 <- dat[,cor.test(Beta, MR.OR)]
    data.table(PRS = this_prs, endpoint = this_endpoint, r = c1$estimate,
               L95 = c1$conf.int[1], U95 = c1$conf.int[2], P = c1$p.value)
  }
}
cor_assocs[, label := sprintf("Pearson r = %s\n[95%% CI: %s-%s]\nP=%s",
             round(r*100)/100, round(L95*100)/100, round(U95*100)/100,
             my_format_pval(P))]

g <- ggplot(comp[anno != "drop"]) +
  aes(x = MR.OR, xmin = MR.OR.L95, xmax = MR.OR.U95,
      y = HR, ymin = HR.L95, ymax = HR.U95,
      color = anno) +
  geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.8, size=0.5) +
  geom_errorbar(width=0, alpha=0.8, size=0.5) +
  geom_point(shape = 19, size=1.3) +
  facet_wrap(~ PRS, nrow=1) +
  scale_x_log10(name = "Causal OR (95% CI)") +
  scale_y_log10(name = "Protein HR (95% CI)") +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g, width=5.8, height=2.4, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/MR_HES_compare.pdf")




