# Update the Olink analysis
library(data.table)
library(foreach)
library(ggplot2)
library(ggrepel)
library(ggrastr)

out_dir <- "analyses/pub/cardiometabolic_proteins/review1"

# Load in associations and filter to those with FDR < 0.05
soma_assocs <- fread("analyses/pub/cardiometabolic_proteins/all_assocs.tsv")
soma_assocs <- soma_assocs[Prot.FDR < 0.05]

# Load somalogic and olink levels, then split out samples who have 
# both measurements
soma_prot <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")
olink_prot <- fread("analyses/processed_traits/olink_proteins/traits.tsv")
olink_has_soma <- olink_prot[IID %in% unique(soma_prot$IID)]
olink_prot <- olink_prot[!(IID %in% unique(soma_prot$IID))]
soma_has_olink <- soma_prot[IID %in% unique(olink_has_soma$IID)]

# Filter to proteins measured on both platforms:
soma_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv", na.strings=c("NA", ""))
olink_info <- fread("analyses/processed_traits/olink_proteins/trait_info.tsv")
olink_info <- olink_info[panel != "neu"] # agreement not in place to use Neurological panel.

common_prot <- intersect(soma_assocs$UniProt, olink_info$UniProt)

soma_info <- soma_info[UniProt.Id.Current.at.Uniprot %in% common_prot]
olink_info <- olink_info[UniProt %in% common_prot]
soma_assocs <- soma_assocs[UniProt %in% common_prot]
olink_prot <- olink_prot[variable %in% olink_info$variable]
olink_has_soma <- olink_has_soma[variable %in% olink_info$variable]
soma_has_olink <- soma_has_olink[variable %in% soma_info$variable]

# Add protein names (Gene symbols)
olink_prot[olink_info, on = .(variable), UniProt := UniProt]
olink_prot[unique(soma_assocs[,.(UniProt, Gene)]), on = .(UniProt), Gene := i.Gene]
olink_prot <- olink_prot[,.(IID, Gene, value)]

olink_has_soma[olink_info, on = .(variable), UniProt := UniProt]
olink_has_soma[unique(soma_assocs[,.(UniProt, Gene)]), on = .(UniProt), Gene := i.Gene]
olink_has_soma <- olink_has_soma[,.(IID, Gene, value)]

soma_has_olink[soma_info, on = .(variable), Gene := Gene.Name]
soma_has_olink <- soma_has_olink[, .(IID, Gene, variable, value)]

# Compute and output cohort characteristics for the two sets of Olink samples
pheno <- fread("analyses/processed_traits/phenotypes.tsv", na.strings=c("NA", ""))
pheno[wt_bl == 777, wt_bl := NA] # bad coding
pheno[, bmi_raw := wt_bl/ht_bl^2]
pheno[, bmi := bmi_raw]
pheno[ht_bl < 1.47, bmi := NA_real_] # clinical cutoff for dwarfism
pheno[ht_bl > 2.1, bmi := NA_real_] # clinical cutoff for gigantism
pheno[wt_bl < 50 | wt_bl > 160, bmi := NA_real_] # NHS restrictions for weight

olink_pheno <- pheno[(IID %in% unique(olink_prot$IID))]
sink(sprintf("%s/olink_cohort_characteristics.txt", out_dir))
print("Total N:")
olink_pheno[, .N]
print("Age:")
olink_pheno[, summary(agePulse)]
print("Sex (2 = woman):")
olink_pheno[, table(sexPulse)]
print("Height (self-reported, m):")
olink_pheno[!is.na(bmi), summary(ht_bl)]
print("Weight (self-reported, kg):")
olink_pheno[!is.na(bmi), summary(wt_bl)]
print("BMI:")
olink_pheno[, summary(bmi)]
print("Failing BMI qc:")
olink_pheno[is.na(bmi), .N]
sink()

overlap_pheno <- pheno[(IID %in% unique(olink_has_soma$IID))]
sink(sprintf("%s/olink_somalogic_overlap_cohort_characteristics.txt", out_dir))
print("Total N:")
overlap_pheno[, .N]
print("Age:")
overlap_pheno[, summary(agePulse)]
print("Sex (2 = woman):")
overlap_pheno[, table(sexPulse)]
print("Height (self-reported, m):")
overlap_pheno[!is.na(bmi), summary(ht_bl)]
print("Weight (self-reported, kg):")
overlap_pheno[!is.na(bmi), summary(wt_bl)]
print("BMI:")
overlap_pheno[, summary(bmi)]
print("Failing BMI qc:")
overlap_pheno[is.na(bmi), .N]
sink()

# Load PRS levels and genetic PCS
prs <- rbind(idcol = "PRS",
  "Coronary Artery Disease"=fread("analyses/GRS_profiles/CAD_metaGRS/profile.sscore.gz"),
  "Type 2 Diabetes"=fread("analyses/GRS_profiles/T2D_2018/profile.sscore.gz")
)

pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Fit associations with olink proteins in the independent samples
olink_assocs <- foreach(gene_id = unique(olink_prot$Gene), .combine=rbind) %do% {
  foreach(prs_id = soma_assocs[Gene == gene_id, unique(PRS)], .combine=rbind) %do% {
    dat <- prs[PRS == prs_id]
    dat <- dat[pcs, on = .(IID)]
    dat <- dat[olink_prot[Gene == gene_id], on = .(IID)]
    l1 <- lm(value ~ scale(score_sum) + PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10, data=dat)
    ci <- confint(l1)
		data.table(PRS=prs_id, Gene=gene_id, Beta=coef(summary(l1))[2,1],
		  				 L95=ci[2,1], U95=ci[2,2], Pvalue=coef(summary(l1))[2,4])
  }
}

# To compare protein levels, we need to adjust them for additional technical covariates.
batch <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")

soma_has_olink <- soma_has_olink[pcs, on = .(IID), nomatch = 0]
soma_has_olink <- soma_has_olink[batch, on = .(IID), nomatch = 0]
soma_has_olink[, value := lm(value ~ factor(batch) + PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=.(variable, Gene)]
soma_has_olink <- soma_has_olink[, .(value = mean(value)), by = .(IID, Gene)]

olink_has_soma <- olink_has_soma[pcs, on = .(IID), nomatch = 0]
olink_has_soma[, value := lm(value ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=.(Gene)]
olink_has_soma <- olink_has_soma[, .(IID, Gene, value)]

# Build comparison tables
soma_assocs <- soma_assocs[,.(PRS, Gene, Beta=Prot.Beta, L95=Prot.L95, U95=Prot.U95, Pvalue=Prot.Pvalue)]
soma_assocs <- unique(soma_assocs)

assoc_comp <- merge(soma_assocs, olink_assocs, by=c("PRS", "Gene"), suffixes=c(".soma", ".oli"))

prot_comp <- merge(soma_has_olink, olink_has_soma, by=c("IID", "Gene"), suffixes=c(".soma", ".oli"))
prot_comp[, Gene := factor(Gene, levels=assoc_comp[order(Beta.soma), Gene])]

# Plot comparison of associations in the independent cohorts:
plot_min <- assoc_comp[, min(c(L95.soma, L95.oli))]
plot_max <- assoc_comp[, max(c(U95.soma, U95.oli))]
expand <- (plot_max - plot_min) * 0.1
plot_min <- plot_min - expand
plot_max <- plot_max + expand

rdt <- data.table(   x = c(plot_min, 0, 0, 0, plot_max),
                  ymin = c(0, 0, 0, plot_min, plot_min),
                  ymax = c(plot_max, plot_max, 0, 0, 0))

col_map <- c(
  "Coronary Artery Disease"="#7570b3",
  "Type 2 Diabetes"="#d95f02",
  "Chronic Kidney Disease"="#1b9e77"
)

cor_assocs <- assoc_comp[, cor.test(Beta.soma, Beta.oli)]
cor_label <- sprintf("Pearson r = %s\n[95%% CI: %s-%s]", 
                     round(cor_assocs$estimate*100)/100,
                     round(cor_assocs$conf.int[1]*100)/100,
                     round(cor_assocs$conf.int[2]*100)/100)

g1 <- ggplot(assoc_comp) +
  aes(x = Beta.soma, y = Beta.oli,
      xmin = L95.soma, ymin = L95.oli,
      xmax = U95.soma, ymax = U95.oli,
      label = Gene,
      color=PRS, fill=PRS) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
  geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#737373") +
  geom_errorbarh(height=0, alpha=1, size=0.5) +
  geom_errorbar(width=0, alpha=1, size=0.5) +
  geom_point(shape = 21, size=1.3, color="#00000000") +
  geom_text_repel(color="black", size=2, nudge_x=0.01) +
  geom_text(inherit.aes=FALSE, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1, label = cor_label) +
  scale_x_continuous(name = "Beta (95% CI) somalogic platform", expand=c(0,0)) +
  scale_y_continuous(name = "Beta (95% CI) olink platform", expand=c(0,0)) +
  scale_color_manual(name = "PRS", values=col_map,
                     guide = guide_legend(direction="horizontal", nrow=3)) +
  scale_fill_manual(name = "PRS", values=col_map,
                    guide = guide_legend(direction="horizontal", nrow=3)) +
  coord_fixed(ratio=1) +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g1, width=3.6, height=3.6, file=sprintf("%s/soma_olink_independent_association_compare.pdf", out_dir), useDingbats=FALSE)

# Compare protein levels in both panels
cor_label <- prot_comp[, .(
  label = sprintf("Pearson r = %s\n[95%% CI: %s-%s]",
                  round(cor.test(value.soma, value.oli)$estimate*100)/100,
                  round(cor.test(value.soma, value.oli)$conf.int[1]*100)/100,
                  round(cor.test(value.soma, value.oli)$conf.int[2]*100)/100)
), by=Gene]

g2 <- ggplot(prot_comp, aes(x=value.soma, y=value.oli)) +
  geom_hline(yintercept=0, size=0.3, colour="#bdbdbd") +
  geom_vline(xintercept=0, size=0.3, colour="#bdbdbd") +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#bdbdbd", size=0.3) +
  geom_point_rast(alpha=0.3, raster.width=3, raster.height=3, size=1) +
  geom_text(data=cor_label, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1, aes(label=label)) +
  facet_wrap(~ Gene, ncol=2) + 
  scale_x_continuous(name="SOMAscan levels", breaks=c(-3, 0, 3), limits=c(-4, 4)) +
  scale_y_continuous(name="Olink levels", breaks=c(-3, 0, 3), limit=c(-4,4)) +
  theme_bw() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),
        strip.text=element_text(size=10, face="bold"),
        strip.background=element_blank(), panel.grid=element_blank())
ggsave(g2, width=3.6, height=4, file=sprintf("%s/soma_olink_prot_compare.pdf", out_dir))
