library(data.table)
library(foreach)
library(ggplot2)
source("src/utilities/format_pval.R")

# Load PRS to protein associations and protein to disease associations
prs_prot <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prot_hes <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")

# Extract protein level associations
prs_prot <- unique(prs_prot[,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])
prot_hes <- unique(prot_hes[,.(endpoint, Target, UniProt, Gene, HR, HR.L95, HR.U95, HR.P, HR.FDR)])

# Map between outcomes
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  endpoint = c("Atrial fibrillation", "Myocardial infarction", "End stage renal disease", "Ischaemic stroke", "Diabetes")
)

prs_prot[map, on = .(PRS), endpoint := endpoint]
prot_hes[map, on = .(endpoint), PRS := PRS]

# Filter to proteins associated with PRS at FDR < 0.05 and endpoints with >= 10 events.
prs_prot <- prs_prot[PRS.FDR < 0.05]
prot_hes <- prot_hes[endpoint %in% c("Atrial fibrillation", "Myocardial infarction", "Diabetes")]

# Build comparison table
comp <- merge(prs_prot, prot_hes, by=c("PRS", "endpoint", "Target", "UniProt", "Gene"))

# Get correlation labels
tests <- unique(comp[,.(PRS, endpoint)])
cor_assocs <- foreach(test_idx = tests[,.I], .combine=rbind) %do% {
  this_prs <- tests[test_idx, PRS]
  this_end <- tests[test_idx, endpoint]

  dat <- comp[PRS == this_prs & endpoint == this_end]
  if (nrow(dat) > 1) {
    c1 <- dat[,cor.test(PRS.Beta, log(HR))]
    data.table(PRS = this_prs, endpoint = this_end, r = c1$estimate,
               L95 = c1$conf.int[1], U95 = c1$conf.int[2], P = c1$p.value)
  }
}
cor_assocs[, label := sprintf("Pearson r = %s\n[95%% CI: %s-%s]\nP=%s",
             round(r*100)/100, round(L95*100)/100, round(U95*100)/100,
             my_format_pval(P))]

# Plot comparison
comp[, sig_hr := HR.P < 0.05]
comp[, sign := sign(PRS.Beta)]

g <- ggplot(comp) +
  aes(x = PRS.Beta, xmin = PRS.L95, xmax = PRS.U95, 
      y = HR, ymin = HR.L95, ymax = HR.U95,
      shape = sig_hr, color = factor(sign)) +
  geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.5, size=0.5) +
  geom_errorbar(width=0, alpha=0.5, size=0.5) +
  geom_point(size = 1.3) +
  scale_shape_manual(guide=FALSE, values=c("TRUE"=19, "FALSE"=1)) +
  scale_color_manual(guide=FALSE, values=c("-1"="#313695", "1"="#a50026")) +
  geom_text(data=cor_assocs, inherit.aes=FALSE, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1, aes(label=label)) +
  facet_wrap(~ PRS) + 
  scale_x_continuous(name = "PRS Beta (95% CI)", breaks=c(-.1, 0, .1)) +
  scale_y_log10(name = "Protein HR (95% CI)", breaks=c(0.3, 0.5, 1, 2, 3)) +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g, width=7.2/2, height=2, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/prs_prot_hes_compare.pdf")


