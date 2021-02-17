library(data.table)
library(foreach)
library(doMC)
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(lubridate)

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
prev_24m <- unique(hes[prevalent_24m == 1 & phenotype %in% cardiometabolic, .(IID)]) 

# Load PRS
prs <- rbind(idcol="PRS",
  CAD_PRS = fread("analyses/GRS_profiles/CAD_metaGRS/profile.sscore.gz"),
  IS_PRS = fread("analyses/GRS_profiles/Stroke_metaGRS/profile.sscore.gz"),
  T2D_PRS = fread("analyses/GRS_profiles/T2D_2018/profile.sscore.gz"),
  AF_PRS = fread("analyses/GRS_profiles/Afib_2018/profile.sscore.gz"),
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

# Load phenotype data
pheno <- fread("analyses/processed_traits/phenotypes.tsv")
pheno <- pheno[!is.na(IID)]

# Compute age at follow-up
pheno <- pheno[,.(IID, sex = sexPulse, age_bl = agePulse, ht_bl, wt_bl, attendanceDate, attendanceDate_24m)]
pheno[, attendanceDate := as.IDate(attendanceDate, format="%d%B%Y")]
pheno[, attendanceDate_24m := as.IDate(attendanceDate_24m, format="%d%B%Y")]
pheno[, follow_24m := time_length(as.Date(attendanceDate_24m) - as.Date(attendanceDate), unit="year")]
pheno[, age_24m := age_bl + round(follow_24m, digits=1)]

# Compute BMI (only available for baseline)
pheno[wt_bl == 777, wt_bl := NA] # bad coding
pheno[, bmi_bl := wt_bl/ht_bl^2]
pheno[ht_bl < 1.47, bmi_bl := NA_real_] # clinical cutoff for dwarfism
pheno[ht_bl > 2.1, bmi_bl := NA_real_] # clinical cutoff for gigantism
pheno[wt_bl < 50 | wt_bl > 160, bmi_bl := NA_real_] # NHS restrictions for weight

# Filter columns
pheno <- pheno[,.(IID, sex, age_bl, bmi_bl, age_24m)]

# Flag prevalent cases
pheno[, prev_bl := FALSE]
pheno[prev, on = .(IID), prev_bl := TRUE]
pheno[, prev_24m := ifelse(is.na(age_24m), NA, FALSE)]
pheno[prev_24m, on = .(IID), prev_24m := TRUE]

# Load SomaLogic aptamer levels
soma <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")

# Adjust for batch
batch <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")
soma <- soma[batch, on = .(IID), nomatch=0]
soma[, value := lm(value ~ factor(batch))$residuals, by=variable]
soma[, soma_ivt_adj_batch := scale(value)]
soma <- soma[,.(IID, variable, soma_ivt_adj_batch)]

# Drop prevalent cases
prs <- prs[!prev, on = .(IID)]
soma <- soma[!prev, on = .(IID)]

# Load protein information.
soma_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

# Filter to aptamers passing QC
soma_info <- soma_info[Type == "Protein"]

# Select columns
soma_info <- soma_info[Type == "Protein", .(variable, SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot,
                                    Gene=Gene.Name, Entrez_id=Entrez.Gene.ID, chr, start, end)]

# Fix bad entries (Aptamers for the same target with different/missing gene/uniprot information)
soma_info[Target == "14-3-3 protein family", UniProt := "P61981|Q04917"]
soma_info[Target == "Induced myeloid leukemia cell differentiation protein Mcl-1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("Q07820", "MCL1", "4170", "1", "150547027", "150552214")]
soma_info[Target == "Protein delta homolog 1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("P80370", "DLK1", "8788", "14", "101193202", "101201467")]
soma_info[Target == "Stromal cell-derived factor 1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("P48061", "CXCL12", "6387", "10", "44865601", "44880545")]

# Filter to aptamers passing QC
soma <- soma[variable %in% soma_info$variable]

# Load in PRS to protein associations and filter to those with FDR < 0.05
soma_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
soma_assocs <- soma_assocs[FDR < 0.05]

# Load in olink qPCR data
olink_qpcr <- fread("analyses/processed_traits/olink_proteins/traits.tsv")

# Remove people with prevalent cardiometabolic disease at the olink timepoint
olink_qpcr <- olink_qpcr[!prev_24m, on = .(IID)]

# Flag samples independent of the somalogic measures and samples with somalogic measures
olink_qpcr[, has_soma := FALSE]
olink_qpcr[IID %in% unique(soma$IID), has_soma := TRUE]

# Load olink information sheet and filter to proteins also measured on the somalogic platform
olink_qpcr_info <- fread("analyses/processed_traits/olink_proteins/trait_info.tsv")
olink_qpcr_info <- olink_qpcr_info[panel != "neu"] # agreement not in place to use Neurological panel.
common_prot <- intersect(unique(soma_assocs$UniProt), olink_qpcr_info$UniProt)

# Compute PRS to olink protein associations for these proteins in both sub-cohorts
common_tests <- unique(soma_assocs[UniProt %in% common_prot, .(PRS, UniProt)])
olink_qpcr_assocs <- foreach(somaset = c(TRUE, FALSE), .combine=rbind) %:%
  foreach(tidx = common_tests[,.I], .combine=rbind) %do% {
    this_prs <- common_tests[tidx, PRS]
    this_up <- common_tests[tidx, UniProt]
    this_gene <- soma_info[UniProt == this_up, unique(Gene)]
    this_target <- soma_info[UniProt == this_up, unique(Target)]
    this_ovar <- olink_qpcr_info[UniProt == this_up, variable]

    dat <- olink_qpcr[has_soma == somaset & variable == this_ovar]
    dat <- merge(dat, prs[PRS == this_prs], by="IID")

    l1 <- lm(scale(value) ~ scale(prs_adj_pcs), data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)
    data.table(has_soma = somaset, PRS = this_prs, Target = this_target, UniProt = this_up, Gene=this_gene, 
               Beta = cf[2,1], SE = cf[2,2], L95 = ci[2,1], U95 = ci[2,2], P = cf[2,4])
}

# Plot comparison of PRS to protein associations in both panels
soma_assocs <- unique(soma_assocs[,.(PRS, Target, UniProt, Gene, Beta, SE, L95, U95, P, FDR)])
comp_qpcr_assocs <- merge(soma_assocs, olink_qpcr_assocs, by= c("PRS", "Target", "UniProt", "Gene"), suffixes=c("", ".olink"))

fwrite(comp_qpcr_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/olink_qpcr_prs_assoc_compare.txt")

col_map <- c(
	"CAD_PRS"="#7570b3",
	"T2D_PRS"="#d95f02",
	"CKD_PRS"="#1b9e77"
)

for (somaset in c(TRUE, FALSE)) {
	assoc_comp <- comp_qpcr_assocs[has_soma == somaset]

	plot_min <- assoc_comp[, min(c(L95, L95.olink))]
	plot_max <- assoc_comp[, max(c(U95, U95.olink))]
	expand <- (plot_max - plot_min) * 0.1
	plot_min <- plot_min - expand
	plot_max <- plot_max + expand

	rdt <- data.table(   x = c(plot_min, 0, 0, 0, plot_max),
										ymin = c(0, 0, 0, plot_min, plot_min),
										ymax = c(plot_max, plot_max, 0, 0, 0))

	cor_assocs <- assoc_comp[, cor.test(Beta, Beta.olink)]
	cor_label <- sprintf("Pearson r = %s\n[95%% CI: %s-%s]",
											 round(cor_assocs$estimate*100)/100,
											 round(cor_assocs$conf.int[1]*100)/100,
											 round(cor_assocs$conf.int[2]*100)/100)

	g <- ggplot(assoc_comp) +
		aes(x = Beta, y = Beta.olink,
				xmin = L95, ymin = L95.olink,
				xmax = U95, ymax = U95.olink,
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
		scale_x_continuous(name = "Beta (95% CI)logic platform", expand=c(0,0)) +
		scale_y_continuous(name = "Beta (95% CI).olinknk platform", expand=c(0,0)) +
		scale_color_manual(name = "PRS", values=col_map,
											 guide = guide_legend(direction="horizontal", nrow=3)) +
		scale_fill_manual(name = "PRS", values=col_map,
											guide = guide_legend(direction="horizontal", nrow=3)) +
		coord_fixed(ratio=1) +
		theme_bw() + theme(
			axis.title=element_text(size=8), axis.text=element_text(size=8),
			panel.grid=element_blank(), legend.position="bottom"
		)

   if (somaset) {
     fname <- "analyses/pub/cardiometabolic_proteins/review2/soma_olink_overlap_association_compare.pdf"
   } else { 
     fname <- "analyses/pub/cardiometabolic_proteins/review2/soma_olink_independent_association_compare.pdf"
   }

	ggsave(g, width=3.6, height=3.6, useDingbats=FALSE, file=fname)
}

# Compare protein levels in people with both sets of measurements
olink_qpcr[olink_qpcr_info, on = .(variable), UniProt := UniProt]
soma <- soma[soma_info, on = .(variable), UniProt := UniProt]

olink_qpcr <- olink_qpcr[UniProt %in% common_prot]
soma <- soma[UniProt %in% common_prot]

soma <- soma[, .(value = mean(soma_ivt_adj_batch)), by=.(IID, UniProt)]

comp_prot <- merge(olink_qpcr, soma, by=c("IID", "UniProt"), suffixes=c(".soma", ".olink"))
comp_prot[soma_info, on = .(UniProt), Gene := Gene]

cor_label <- comp_prot[, .(
  label = sprintf("Pearson r = %s\n[95%% CI: %s-%s]",
                  round(cor.test(value.soma, value.olink)$estimate*100)/100,
                  round(cor.test(value.soma, value.olink)$conf.int[1]*100)/100,
                  round(cor.test(value.soma, value.olink)$conf.int[2]*100)/100)
), by=Gene]

g <- ggplot(comp_prot, aes(x=value.soma, y=value.olink)) +
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
ggsave(g, width=3.6, height=4, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/soma_vs_olink_levels.pdf")

# Load Olink NovaSeq data
olink_nova <- fread("analyses/processed_traits/olink_ngs/traits.tsv")

# Drop people with prevalent cardiometabolic disease
olink_nova <- olink_nova[!prev_24m, on = .(IID)]

# Filter to proteins measured on both platforms
olink_nova_info <- fread("analyses/processed_traits/olink_ngs/trait_info.tsv")
up_map <- olink_nova_info[,.(UniProt = strsplit(UniProt, "_")[[1]]), by=.(Panel, variable)]
common_prot <- intersect(soma_assocs$UniProt, up_map$UniProt)

# Test each Olink protein for association with PRS
common_tests <- soma_assocs[UniProt %in% common_prot]
olink_nova_assocs <- foreach(tidx = common_tests[,.I], .combine=rbind) %do% {
	this_prs <- common_tests[tidx, PRS]
	this_up <- common_tests[tidx, UniProt]
	this_gene <- soma_info[UniProt == this_up, unique(Gene)]
	this_target <- soma_info[UniProt == this_up, unique(Target)]
	this_ovar <- olink_nova_info[UniProt == this_up, variable]

	dat <- olink_nova[variable == this_ovar]
	dat <- merge(dat, prs[PRS == this_prs], by="IID")

	l1 <- lm(scale(value) ~ scale(prs_adj_pcs), data=dat)
	cf <- coef(summary(l1))
	ci <- confint(l1)
	data.table(PRS = this_prs, Target = this_target, UniProt = this_up, Gene=this_gene, 
						 Beta = cf[2,1], SE = cf[2,2], L95 = ci[2,1], U95 = ci[2,2], P = cf[2,4])
}

# Plot comparison of associations
assoc_comp <- merge(soma_assocs, olink_nova_assocs, by=c("PRS", "Target", "UniProt", "Gene"), suffixes=c("", ".olink"))
fwrite(assoc_comp, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/olink_ngs_prs_assoc_compare.txt")

plot_min <- assoc_comp[, min(c(L95, L95.olink))]
plot_max <- assoc_comp[, max(c(U95, U95.olink))]
expand <- (plot_max - plot_min) * 0.1
plot_min <- plot_min - expand
plot_max <- plot_max + expand

rdt <- data.table(   x = c(plot_min, 0, 0, 0, plot_max),
									ymin = c(0, 0, 0, plot_min, plot_min),
									ymax = c(plot_max, plot_max, 0, 0, 0))

cor_assocs <- assoc_comp[, cor.test(Beta, Beta.olink)]
cor_label <- sprintf("Pearson r = %s\n[95%% CI: %s-%s]",
										 round(cor_assocs$estimate*100)/100,
										 round(cor_assocs$conf.int[1]*100)/100,
										 round(cor_assocs$conf.int[2]*100)/100)

g <- ggplot(assoc_comp) +
	aes(x = Beta, y = Beta.olink,
			xmin = L95, ymin = L95.olink,
			xmax = U95, ymax = U95.olink,
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
	scale_x_continuous(name = "Beta (95% CI)logic platform", expand=c(0,0)) +
	scale_y_continuous(name = "Beta (95% CI).olinknk platform", expand=c(0,0)) +
	scale_color_manual(name = "PRS", values=col_map,
										 guide = guide_legend(direction="horizontal", nrow=3)) +
	scale_fill_manual(name = "PRS", values=col_map,
										guide = guide_legend(direction="horizontal", nrow=3)) +
	coord_fixed(ratio=1) +
	theme_bw() + theme(
		axis.title=element_text(size=8), axis.text=element_text(size=8),
		panel.grid=element_blank(), legend.position="bottom"
	)

ggsave(g, width=3.6, height=3.6, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/olink_nova_soma_assoc_compare.pdf")

# Compare Olink NovaSeq to Olink qPCR proteins
common_prot <- unique(olink_qpcr_assocs$UniProt)
olink_nova[olink_nova_info, on = .(variable), UniProt := i.UniProt]
olink_nova <- olink_nova[UniProt %in% common_prot]

comp_prot <- merge(olink_qpcr, olink_nova, by=c("IID", "UniProt"), suffixes=c(".qpcr", ".nova"))
comp_prot[soma_info, on = .(UniProt), Gene := Gene]

cor_label <- comp_prot[, .(
  label = sprintf("Pearson r = %s\n[95%% CI: %s-%s]",
                  round(cor.test(value.qpcr, value.nova)$estimate*100)/100,
                  round(cor.test(value.qpcr, value.nova)$conf.int[1]*100)/100,
                  round(cor.test(value.qpcr, value.nova)$conf.int[2]*100)/100)
), by=Gene]

g <- ggplot(comp_prot, aes(x=value.qpcr, y=value.nova)) +
  geom_hline(yintercept=0, size=0.3, colour="#bdbdbd") +
  geom_vline(xintercept=0, size=0.3, colour="#bdbdbd") +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#bdbdbd", size=0.3) +
  geom_point_rast(alpha=0.3, raster.width=3, raster.height=3, size=1) +
  geom_text(data=cor_label, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1, aes(label=label)) +
  facet_wrap(~ Gene, ncol=2) +
  scale_x_continuous(name="Olink qPCR levels", breaks=c(-3, 0, 3), limits=c(-4, 4)) +
  scale_y_continuous(name="Olink NGS levels", breaks=c(-3, 0, 3), limit=c(-4,4)) +
  theme_bw() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),
        strip.text=element_text(size=10, face="bold"),
        strip.background=element_blank(), panel.grid=element_blank())
ggsave(g, width=3.6, height=4, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/olink_qpcr_vs_nova_levels.pdf")




