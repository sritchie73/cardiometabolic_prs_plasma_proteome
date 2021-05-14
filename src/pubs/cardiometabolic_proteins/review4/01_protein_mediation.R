library(data.table)
library(ggplot2)

# Load mediation analysis results
med <- fread("analyses/pub/cardiometabolic_proteins/review3/protein_mediation.txt")

# split into effect types
nde <- med[Mediation == "natural direct effect"] 
nie <- med[Mediation == "natural indirect effect"] 
te <- med[Mediation == "total effect"]

# Compute % mediated and 95% CIs
nde[te, on = .(PRS, Target, UniProt, Gene, Endpoint, Aptamer), PTE.L95 := L95 / i.logOR]
nde[te, on = .(PRS, Target, UniProt, Gene, Endpoint, Aptamer), PTE.U95 := U95 / i.logOR]
nde[te, on = .(PRS, Target, UniProt, Gene, Endpoint, Aptamer), PTE.L95.aptamer := L95.aptamer / i.logOR.aptamer]
nde[te, on = .(PRS, Target, UniProt, Gene, Endpoint, Aptamer), PTE.U95.aptamer := U95.aptamer / i.logOR.aptamer]
nde <- nde[, .(PRS, Target, UniProt, Gene, Endpoint, logOR, SE, L95, U95, PTE, PTE.L95, PTE.U95, P,
               Aptamer, logOR.aptamer, SE.aptamer, L95.aptamer, U95.aptamer, PTE.aptamer,
               PTE.L95.aptamer, PTE.U95.aptamer, P.aptamer)]

nie[te, on = .(PRS, Target, UniProt, Gene, Endpoint, Aptamer), PTE.L95 := L95 / i.logOR]
nie[te, on = .(PRS, Target, UniProt, Gene, Endpoint, Aptamer), PTE.U95 := U95 / i.logOR]
nie[te, on = .(PRS, Target, UniProt, Gene, Endpoint, Aptamer), PTE.L95.aptamer := L95.aptamer / i.logOR.aptamer]
nie[te, on = .(PRS, Target, UniProt, Gene, Endpoint, Aptamer), PTE.U95.aptamer := U95.aptamer / i.logOR.aptamer]
nie <- nie[, .(PRS, Target, UniProt, Gene, Endpoint, logOR, SE, L95, U95, PTE, PTE.L95, PTE.U95, P,
               Aptamer, logOR.aptamer, SE.aptamer, L95.aptamer, U95.aptamer, PTE.aptamer,
               PTE.L95.aptamer, PTE.U95.aptamer, P.aptamer)]

te <- te[, .(PRS, Target, UniProt, Gene, Endpoint, logOR, SE, L95, U95, PTE=1, PTE.L95=1, PTE.U95=1, P,
             Aptamer, logOR.aptamer, SE.aptamer, L95.aptamer, U95.aptamer, PTE.aptamer=1,
             PTE.L95.aptamer=1, PTE.U95.aptamer=1, P.aptamer)]

# Combine and write out
med <- rbind(nde, nie, te)
fwrite(med, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/protein_mediation.txt")

# Redo Fig 2d.
nie[, sig := P < 0.05]
nie <- nie[order(-logOR)][order(-sig)][order(-PRS)]
nie[Gene == "SHBG" & PRS == "CAD_PRS", Gene := "SHBG."]
nie[, Gene := factor(Gene, levels=unique(Gene))]
nie <- unique(nie[,.(PRS, Gene, PTE, PTE.L95, PTE.U95, sig)])

g <- ggplot(nie) +
  aes(x = Gene, y = PTE, ymin = PTE.L95, ymax = PTE.U95, fill=sig) +
  geom_hline(yintercept=0, linetype=2, size=0.2) +
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
ggsave(g, width=110, height=44, units="mm", file="analyses/pub/cardiometabolic_proteins/review4/causal_or_plot.pdf")

# Write out columns for Extended Data Table 2.
old <- fread("analyses/pub/cardiometabolic_proteins/review3/compare_estimates.txt")

format_pct <- function(x) {
  ifelse(x < 0.1,
      sprintf("%.01f%%", round(x * 1000)/10),
      sprintf("%d%%", round(x * 100)))
}

format_ci <- function(l95, u95) {
  sprintf("[%s, %s]", format_pct(l95), format_pct(u95))
}

nie[Gene == "SHBG.", Gene := "SHBG"]
old[nie, on = .(PRS, Gene), Causal.PTE.95CI := format_ci(PTE.L95, PTE.U95)]

new <- old[, .(PRS, Endpoint, Gene, PRS.Beta, PRS.95CI, PRS.P, PRS.FDR, Polygenicity,
               HR, HR.95CI, HR.P, Causal.PTE, Causal.PTE.95CI, Causal.P)]
fwrite(new, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/compare_estimates.txt")
