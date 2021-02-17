library(data.table)
library(ggplot2)
library(cowplot)

system("mkdir -p analyses/pub/cardiometabolic_proteins/review3", wait=TRUE)

# Load summary statistics generated under review2/
prs_hrs <- fread("analyses/pub/cardiometabolic_proteins/review2/prs_hes_associations.txt")
ckd_egfr <- fread("analyses/pub/cardiometabolic_proteins/review2/ckd_prs_egfr_associations.txt")

# Filter to analyses done in the sub-cohort with somalogic data
prs_hrs <- prs_hrs[Cohort == "Soma"]
ckd_egfr <- ckd_egfr[Cohort == "Soma"]

# Drop models with < 10 events
prs_hrs <- prs_hrs[events >= 10]

# Create plots
g1 <- ggplot(prs_hrs) +
  aes(x=HR, xmin=HR.L95, xmax=HR.U95, y=factor(test)) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0, alpha=0.8) +
  geom_point(shape=18, size=3) +
  xlab("Hazard Ratio per SD increase of PRS") +
  ylab("") +
  theme_bw()

g2 <- ggplot(ckd_egfr) +
  aes(x=beta, xmin=L95, xmax=U95, y=factor(test)) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, alpha=0.8) +
  geom_point(shape=18, size=3) +
  xlab("CKD-EPI eGFR (mL/min/1.73 m2) per SD of PRS") +
  ylab("") +
  theme_bw()

g <- plot_grid(g1, g2, nrow=2, rel_heights=c(0.7, 0.3))
ggsave(g, width=6.2, height=3.5, file="analyses/pub/cardiometabolic_proteins/review3/prs_validation.pdf")


