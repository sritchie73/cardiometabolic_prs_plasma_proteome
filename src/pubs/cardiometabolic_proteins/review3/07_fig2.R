library(data.table)
library(ggplot2)

# Load data
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
hes_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")
med <- fread("analyses/pub/cardiometabolic_proteins/review3/protein_mediation.txt")

# Map between tables
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  endpoint = c("Atrial fibrillation", "Myocardial infarction", "End stage renal disease", "Ischaemic stroke", "Diabetes")
)

prs_assocs[map, on = .(PRS), Endpoint := i.endpoint]
hes_assocs[map, on = .(endpoint), PRS := i.PRS]
setnames(hes_assocs, "endpoint", "Endpoint")

# Filter to significant testable associations
prs_assocs <- unique(prs_assocs[FDR < 0.05,.(PRS, Endpoint, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])
prs_assocs <- prs_assocs[!(PRS %in% c("CKD_PRS", "IS_PRS"))]
hes_assocs <- hes_assocs[prs_assocs[,.(PRS, Endpoint, Gene)], on = .(PRS, Endpoint, Gene)]
hes_assocs <- unique(hes_assocs[, .(PRS, Endpoint, Gene, HR, HR.L95, HR.U95, HR.P)])
med <- unique(med[Mediation == "natural indirect effect", .(PRS, Endpoint, Gene, OR = exp(logOR), L95 = exp(L95), U95 = exp(U95), P, PTE)])

# Plot HR
hes_assocs <- hes_assocs[order(-HR)][order(-PRS)]
hes_assocs[Gene == "SHBG" & PRS == "CAD_PRS", Gene := "SHBG."]
hes_assocs[, Gene := factor(Gene, levels=unique(Gene))]
hes_assocs[, sig := HR.P < 0.05]

g <- ggplot(hes_assocs) +
  aes(x = Gene, y = HR, ymin = HR.L95, ymax = HR.U95, fill=sig) +
  geom_hline(yintercept=1, linetype=2, size=0.2) +
  geom_errorbar(width=0, size=0.2) +
  geom_point(shape=23, stroke=0.4, size=0.8) + 
  scale_fill_manual(guide=FALSE, values=c("TRUE"="black", "FALSE"="white")) +
  facet_grid(~ PRS, scales="free_x", space="free_x") +
  scale_y_log10(name="HR (95% CI)", breaks=c(0.3, 0.5, 1, 2, 3)) +
  xlab("") +
  theme_bw() +
  theme(axis.title=element_text(size=7), axis.text.y=element_text(size=5), axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
        axis.ticks.y=element_line(size=0.3), axis.ticks.x=element_blank(), panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(size=0.3), panel.grid.minor.y=element_blank(),
        strip.text=element_text(size=7), strip.background=element_blank(), panel.spacing=unit(0.5, "mm"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"))
ggsave(g, width=110, height=44, units="mm", file="analyses/pub/cardiometabolic_proteins/review3/hr_plot.pdf")

# Plot causal OR
med[, sig := P < 0.05]
med <- med[order(-OR)][order(-sig)][order(-PRS)]
med[Gene == "SHBG" & PRS == "CAD_PRS", Gene := "SHBG."]
med[, Gene := factor(Gene, levels=unique(Gene))]

g <- ggplot(med) +
  aes(x = Gene, y = OR, ymin = L95, ymax = U95, fill=sig) +
  geom_hline(yintercept=1, linetype=2, size=0.2) +
  geom_errorbar(width=0, size=0.2) +
  geom_point(shape=23, stroke=0.4, size=0.8) + 
  scale_fill_manual(guide=FALSE, values=c("TRUE"="black", "FALSE"="white")) +
  facet_grid(~ PRS, scales="free_x", space="free_x") +
  ylab("Causal OR (95% CI)") + xlab("") +
  theme_bw() +
  theme(axis.title=element_text(size=7), axis.text.y=element_text(size=5), axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
        axis.ticks.y=element_line(size=0.3), axis.ticks.x=element_blank(), panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(size=0.3), panel.grid.minor.y=element_blank(),
        strip.text=element_text(size=7), strip.background=element_blank(), panel.spacing=unit(0.5, "mm"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"))
ggsave(g, width=110, height=44, units="mm", file="analyses/pub/cardiometabolic_proteins/review3/causal_or_plot.pdf")

