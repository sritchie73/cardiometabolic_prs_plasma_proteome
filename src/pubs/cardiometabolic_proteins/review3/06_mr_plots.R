library(data.table)
library(ggplot2)
library(foreach)
source("src/utilities/prot_pval.R")
source("src/07_job_scripts/07_helpers/mr_functions.R")

# Load MR Estimates and filter to PRS-associated proteins
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
mr_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/collated_mr.txt")

prs_assocs <- unique(prs_assocs[,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])
mr_assocs <- unique(mr_assocs[Method == "(Intercept)",.(GWAS, Target, UniProt, Gene, MR=exp(mr_estimate), MR.L95=exp(mr_L95), MR.U95=exp(mr_U95), MR.P=mr_pval, 
                              Pleiotropy=exp(mr_estimate.protein), Pleiotropy.L95=exp(mr_L95.protein), Pleiotropy.U95=exp(mr_U95.protein), Pleiotropy.P=pleiotropy_pval)])

# Map between outcomes
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  GWAS = c("Afib", "CAD", "CKD", "StrokeIS", "T2DadjBMI")
)

prs_assocs[map, on = .(PRS), c("GWAS") := .(GWAS)]
mr_assocs[map, on = .(GWAS), c("PRS") := .(PRS)]

# Filter to FDR < 0.05 PRS to protein associations
prs_assocs <- prs_assocs[PRS.FDR < 0.05]
mr_assocs <- mr_assocs[prs_assocs[,.(PRS, GWAS, Gene)], on = .(PRS, GWAS, Gene), nomatch=0]

# Load instruments
mr_ivs <- fread("analyses/pub/cardiometabolic_proteins/review2/cis_ivs.txt")
mr_ivs[map, on = .(GWAS), PRS := i.PRS]
mr_ivs <- mr_ivs[mr_assocs[,.(PRS, GWAS, Gene)], on = .(PRS, GWAS, Gene)]

# Average effects across aptamers (for IVs used in multiple aptamers)
mr_ivs <- mr_ivs[, .(pQTL.beta=mean(pQTL.beta), pQTL.SE=mean(pQTL.SE), pQTL.P=prot_pvalue(pQTL.P, pQTL.beta)),
                 by=.(PRS, GWAS, Gene, pQTL.chr, pQTL.pos, pQTL.EA, pQTL.OA, gwas.beta, gwas.SE, gwas.P)]

# Create plot label text
mr_assocs <- mr_assocs[order(MR.P)]
mr_assocs[, labtext := sprintf("%s ~ %s", Gene, GWAS)]
mr_assocs[, labtext := factor(labtext, levels=unique(labtext))]
mr_ivs[mr_assocs, on = .(GWAS, Gene), labtext := i.labtext]

# Determine the plot limits
plotlim <- mr_ivs[, .(xmin = min(pQTL.beta - pQTL.SE), xmax = max(pQTL.beta + pQTL.SE),
                      ymin = min(gwas.beta - gwas.SE), ymax = max(gwas.beta + gwas.SE))]

# Make sure to include intercepts
plotlim[xmax < 0, xmax := 0]
plotlim[xmin > 0, xmin := 0]
plotlim[ymax < 0, ymax := 0]
plotlim[ymin > 0, ymin := 0]

# Expand by 5% as per ggplot defaults
expand <- plotlim[, .(xexpand = abs((xmax - xmin))*0.05, yexpand = abs((ymax - ymin))*0.05)]
plotlim[, c("xmin", "xmax", "ymin", "ymax") := .(xmin - expand$xexpand, xmax + expand$xexpand, ymin - expand$yexpand, ymax + expand$yexpand)]

# Now build polygons for 95% CIs. This table is a series of x and y-coordinates for drawing
# the polygon for each protein to disease pair. We need to work out where the polygon intersects 
# the edges of the plot so we know where to define the polygon edges.
mr.ci95 <- foreach(testIdx = mr_assocs[,.I], .combine=rbind) %do% {
	mrtest <- mr_assocs[testIdx]
	intercept.l95 <- mrtest[, log(Pleiotropy.L95)]
	intercept.u95 <- mrtest[, log(Pleiotropy.U95)]
	slope.l95 <- mrtest[, log(MR.L95)]
	slope.u95 <- mrtest[, log(MR.U95)]

	xlim <- plotlim[, c(xmin, xmax)]
	ylim <- plotlim[, c(ymin, ymax)]

	if (intercept.l95 < ylim[1]) intercept.l95 <- ylim[1]
	if (intercept.u95 > ylim[2]) intercept.u95 <- ylim[2]

	poly <- line_ci95_poly(intercept.l95, intercept.u95, slope.l95, slope.u95, xlim, ylim)
	cbind(mr_assocs[testIdx,.(labtext)], poly)
}


# plot
g <- ggplot(mr_ivs) + 
  aes(x=pQTL.beta, xmin=pQTL.beta - pQTL.SE, xmax=pQTL.beta + pQTL.SE,
      y=exp(gwas.beta), ymin=exp(gwas.beta - gwas.SE), ymax=exp(gwas.beta + gwas.SE)) +
  geom_polygon(data=mr.ci95, inherit.aes=FALSE, aes(x=x, y=exp(y)), fill="#fed976", color=NA, size=0) +
  geom_hline(yintercept=1, linetype=3, color="#bdbdbd", size=0.2) + 
  geom_vline(xintercept=0, linetype=3, color="#bdbdbd", size=0.2) + 
  geom_abline(data=mr_assocs, aes(intercept=Pleiotropy, slope=-(1 - MR)), color="#e31a1c", linetype=2, size=0.3) +
  geom_errorbarh(height=0, size=0.2) + geom_errorbar(width=0, size=0.2) +
  geom_point(shape=23, size=0.6, fill="white", stroke=0.4) +
  facet_wrap(~ labtext, nrow=3, dir = "h") + 
  scale_x_continuous(name="Effect of pQTL on protein", expand=c(0,0)) +
  scale_y_continuous(name="OR for pQTL on disease", expand=c(0,0)) + 
  theme_bw() +
  theme(axis.title=element_text(size=7), axis.text=element_text(size=5), panel.grid=element_blank(),
        strip.text=element_blank(), strip.background=element_blank(), panel.spacing=unit(0.5, "mm"),
        axis.ticks=element_line(size=0.3), plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"))
ggsave(g, width=109, height=46, units="mm", file="analyses/pub/cardiometabolic_proteins/review3/mr_dose_response.pdf")


