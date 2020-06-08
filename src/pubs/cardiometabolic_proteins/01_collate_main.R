# Generate tables and figures for publication
library(data.table)
library(foreach)
library(openxlsx)
library(RColorBrewer)
library(pheatmap)
library(seriation)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(ggrastr)
library(gridExtra)
library(NetRep)

source("src/utilities/format_pval.R")
source("src/utilities/prot_pval.R")
source("src/07_job_scripts/07_helpers/mr_functions.R") # for dose response curves

# Explanation for useDingbats=FALSE in ggsave():
# 
# On my Windows machine, the points from scatterplots were being imported
# into inkscape at half their size (and with a small offset, not noticeable unless 
# points are large). This is because R uses Dingbats symbols to render points 
# (see ?pdf) which are then imported incorrectly into inkscape (the point size is 
# 50% of that seen in Adobe when opening the PDF, and points have a slight left offset).
# See also: https://stackoverflow.com/questions/23524262/why-doesnt-inkscape-correctly-read-pdf-files-generated-by-r
# This does not appear to have been an issue in the past when I was using OSX.

# Setup
view_file <- "views/cardiometabolic.txt"
grs_info <- fread(view_file)
GRSs <- grs_info[, GRS_name]
out_dir <- sprintf("analyses/pub/cardiometabolic_proteins")
dir.create(out_dir, recursive = TRUE, showWarnings=FALSE) 

# Load in protein data 
prot <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv", colClasses=c("IID"="character"))
prot_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv", na.strings=c("NA", ""))

# Load in each GRS
score_files <- sprintf("analyses/GRS_profiles/%s/profile.sscore.gz", GRSs)
names(score_files) <- GRSs
levels <- lapply(score_files, fread, colClasses=c("IID"="character"))
levels <- rbindlist(levels, idcol="grs", fill=TRUE)

# -----------------------------------------------------------
# Output cohort summary statistics
# -----------------------------------------------------------
pheno <- fread("analyses/processed_traits/phenotypes.tsv", na.strings=c("NA", ""), colClasses=c("IID"="character"))
pheno[wt_bl == 777, wt_bl := NA] # bad coding
pheno[, bmi_raw := wt_bl/ht_bl^2]
pheno[, bmi := bmi_raw]
pheno[ht_bl < 1.47, bmi := NA_real_] # clinical cutoff for dwarfism
pheno[ht_bl > 2.1, bmi := NA_real_] # clinical cutoff for gigantism
pheno[wt_bl < 50 | wt_bl > 160, bmi := NA_real_] # NHS restrictions for weight

# Sub cohort with protein measurements
prot_pheno <- pheno[IID %in% unique(prot$IID) & IID %in% unique(levels$IID)]

# Output cohort characteristics:
sink(sprintf("%s/cohort_characteristics.txt", out_dir))
print("Total N:")
prot_pheno[, .N]
print("Age:")
prot_pheno[, summary(agePulse)]
print("Sex (2 = woman):")
prot_pheno[, table(sexPulse)]
print("Height (self-reported, m):")
prot_pheno[!is.na(bmi), summary(ht_bl)]
print("Weight (self-reported, kg):")
prot_pheno[!is.na(bmi), summary(wt_bl)]
print("BMI:")
prot_pheno[, summary(bmi)]
print("Failing BMI qc:")
prot_pheno[is.na(bmi), .N]
sink()

# -----------------------------------------------------------
# Plot heatmap of GRS level correlations
# -----------------------------------------------------------

# Load PCs so we can adjust when calculating correlations:
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt", colClasses=c("ID"="character"))
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Filter levels to IIDs with protein measurements, adjust for 10 PCs and convert to wide format
levels_cor <- copy(levels)
levels_cor <- levels_cor[IID %in% unique(prot$IID)]
levels_cor <- levels_cor[pcs, on = .(IID), nomatch=0]
levels_cor <- levels_cor[, .(IID, score = lm(score ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals), by = grs]
levels_cor <- dcast(levels_cor, IID ~ grs, value.var="score")
levels_cor <- as.matrix(levels_cor, rownames=1)

# Calculate pairwise correlation
cor_mat <- cor(levels_cor, method="pearson")
text_mat <- round(cor_mat*100)/100

# Create heatmap
pal <- rev(brewer.pal(name="RdYlBu", n=8))
pal <- c(pal[1:4], "#FFFFFF", pal[5:8])

pheatmap(cor_mat, display_numbers=text_mat,
  color=colorRampPalette(pal)(255), 
  breaks=seq(-1, 1, length=256),
  legend_breaks=seq(-1, 1, length=11),
  fontsize_row=8, fontsize_col=8, fontsize_number=8,
  cellwidth=30, cellheight=30,
  filename=sprintf("%s/grs_cor.pdf", out_dir), width=7, height=6)

# -----------------------------------------------------------
# Collate table of associations
# -----------------------------------------------------------

# Load associations
assocs <- foreach(grs = GRSs) %do% {
  fread(sprintf("analyses/univariate_associations/%s/somalogic_proteins/associations.tsv", grs))
}
names(assocs) <- GRSs
assocs <- rbindlist(assocs, idcol="grs")

# Filter to "good" aptamers
prot_info <- prot_info[Type == "Protein"]

# Add information, filter, average associations across aptamers, and filter at FDR < 0.05:
assocs <- assocs[trait %in% prot_info$variable]
assocs <- assocs[prot_info, on = .(trait=variable), nomatch=0]
prot_assocs <- assocs[, .(beta = mean(beta), l95 = mean(l95), u95 = mean(u95), pval = prot_pvalue(pval, beta)), 
  by = .(grs, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, 
         Gene.Name, Entrez.Gene.ID, chr, start, end, strand)] 
prot_assocs[, fdr := p.adjust(pval, method="fdr"), by=grs]

# Build combined table of relevant information for supp
comb_assocs <- prot_assocs[assocs, 
  on = .(grs, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot,
         Gene.Name, Entrez.Gene.ID, chr, start, end, strand)]
comb_assocs <- comb_assocs[, .(PRS=grs, Target, UniProt, Gene=Gene.Name,
  Entrez.ID = Entrez.Gene.ID, Chr=chr, Start=start, End=end, Strand=strand,
  Prot.Beta=beta, Prot.L95=l95, Prot.U95=u95, Prot.Pvalue=pval, Prot.FDR=fdr, 
  Aptamer=SeqId, Apt.Beta=i.beta, Apt.L95=i.l95, Apt.U95=i.u95, Apt.Pvalue=i.pval, 
  CrossReactivity=Characterization.Info, 
  MassSpecConfirmation=Mass.Spec.Confirmation.in.Matrix)]

# Add in cis-pQTL info
cis_pQTLs <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
cis_pQTLs[prot_info, on = .(SOMAMER_ID), SeqId := SeqId]
comb_assocs[cis_pQTLs, on = .(Aptamer=SeqId), cis_pQTL := "yes"]

# Load in additional mass spec confirmation from Emillson et al.
emilsson <- rbind(fill=TRUE,
  as.data.table(read.xlsx("data/Emilsson_etal_2018/supp_tables.xlsx", sheet="Table S3", startRow=2)),
  as.data.table(read.xlsx("data/Emilsson_etal_2018/supp_tables.xlsx", sheet="Table S4", startRow=2)))
comb_assocs[emilsson, on = .(Gene=Gene.Symbol), MassSpecConfirmation := ifelse(is.na(MassSpecConfirmation), "Emilsson_etal_2018", MassSpecConfirmation)]

# Add pretty GRS name labels
comb_assocs[grs_info, on=.(PRS=GRS_name), PRS := Display_name]

# Organise:
comb_assocs <- comb_assocs[order(Apt.Pvalue)][order(Prot.Pvalue)][order(PRS)]

# Write out 
fwrite(comb_assocs, file=sprintf("%s/all_assocs.tsv", out_dir), sep="\t", quote=FALSE)

# -----------------------------------------------------------
# QQ-plot of all associations
# -----------------------------------------------------------
# Combine info
qq_assocs <- copy(prot_assocs)
qq_assocs[, sig := "FALSE"]
qq_assocs[fdr < 0.1, sig := "Suggestive"]
qq_assocs[fdr < 0.05, sig := "Significant"] 

# Add pretty GRS name labels, using abbrevations because 
# labels too wide otherwise
qq_assocs[grs_info, on=.(grs=GRS_name), Trait := Display_name]

# Generate expected distribution of p-values
qq_assocs <- qq_assocs[order(pval)][order(Trait)]
qq_assocs[, expected := ppoints(.N), by = Trait]

# Transform for plot
qq_assocs[, expected := -log10(expected)]
qq_assocs[, pval := -log10(pval)]

# Get annotations
anno <- qq_assocs[, .SD[1:5], by=grs]
anno[prot_info, on = .(Target=TargetFullName), c("Target", "TargetFullName") := .(i.Target, i.TargetFullName)]
anno <- anno[, .(Trait, Target, TargetFullName, UniProt, Gene.Name, fdr)]

# Plot
g <- ggplot(qq_assocs) +
  aes(x=expected, y=pval) +
  geom_abline(intercept=0, slope=1, linetype=2, color="#737373") +
  geom_line(color="black", size=0.5) +
  geom_point_rast(shape=19, aes(color=sig), size=0.8, raster.width=1.270, raster.height=1.281, raster.dpi=1200) +
  geom_hline(color="#0571b0", linetype=2, yintercept=-log10(0.05)) +
  facet_wrap(~ Trait, nrow=1) +
  scale_x_continuous(name = "Expected -log10 p-values") +
  scale_y_continuous(name = "Observed -log10 p-values") +
  scale_color_manual(values=c("FALSE"="#000000", "Suggestive"="#4575b4", 
                              "Significant"="#d73027"), guide=FALSE) +
  theme_bw() + theme(
    strip.text=element_text(size=8), axis.title=element_text(size=10), 
    axis.text=element_text(size=8), panel.grid=element_blank()
  )
ggsave(g, width=7.2, height=2, file=sprintf("%s/assoc_qq.pdf", out_dir), useDingbats=FALSE)

# -----------------------------------------------------------
# Show heatmap of significant associations 
# -----------------------------------------------------------

# Build data.table for heatmap. Show all proteins where there
# was FDR < 0.05 for any PRS.
heatmap_prots <- prot_assocs[fdr < 0.05, unique(Gene.Name)]
heatmap_dt <- prot_assocs[Gene.Name %in% heatmap_prots]
heatmap_dt[grs_info, on = .(grs=GRS_name), PRS := Display_name]

# Create wide tables for heatmap
beta_mat <- as.matrix(dcast(heatmap_dt, PRS ~ Gene.Name, value.var="beta"), rownames="PRS")
pval_mat <- as.matrix(dcast(heatmap_dt, PRS ~ Gene.Name, value.var="pval"), rownames="PRS")
fdr_mat <- as.matrix(dcast(heatmap_dt, PRS ~ Gene.Name, value.var="fdr"), rownames="PRS")

# Create text matrix
text_mat <- matrix("", nrow=nrow(beta_mat), ncol=ncol(beta_mat),
											dimnames=dimnames(beta_mat))
text_mat[pval_mat < 0.05] <- "o"
text_mat[fdr_mat < 0.1] <- "*"
text_mat[fdr_mat < 0.05] <- "#"

# Order heatmap
prs_order <- c("Chronic Kidney Disease", "Type 2 Diabetes", "Coronary Artery Disease", 
               "Stroke", "Atrial Fibrillation")
prot_order <- foreach(pIdx = seq_along(prs_order), .combine=c) %do% {
  # First show positive the negative associations, keeping aside any proteins
  # associated with the next PRS. Drop any proteins associated with the previous 
  # prs in the order so we don't have duplicates.
  this_prs <- heatmap_dt[PRS == prs_order[pIdx] & fdr < 0.05]
  if (pIdx > 1) {
    assoc_with_last <- heatmap_dt[PRS == prs_order[pIdx - 1] & fdr < 0.05, .(Gene.Name)]
    this_prs <- this_prs[!assoc_with_last, on = .(Gene.Name)]
  }
  if (pIdx != length(prs_order)) {
    assoc_with_next <- heatmap_dt[PRS == prs_order[pIdx + 1] & fdr < 0.05, .(Gene.Name)]
    both_prs <- this_prs[assoc_with_next, on = .(Gene.Name), nomatch = 0]  
    this_prs <- this_prs[!assoc_with_next, on = .(Gene.Name)]
  }
  po <- c(this_prs[beta > 0][order(-beta), Gene.Name], this_prs[beta < 0][order(beta), Gene.Name])
  if (pIdx != length(prs_order)) {
    po <- c(po, both_prs[beta > 0][order(-beta), Gene.Name], both_prs[beta < 0][order(beta), Gene.Name])
  }
  return(po)
}
beta_mat <- beta_mat[prs_order, prot_order]
text_mat <- text_mat[prs_order, prot_order]

# Split heatmap into positive and negative associations
pos_prot <- heatmap_dt[fdr < 0.05 & beta > 0, Gene.Name]
pos_prot <- intersect(prot_order, pos_prot)
neg_prot <- heatmap_dt[fdr < 0.05 & beta < 0, Gene.Name]
neg_prot <- intersect(prot_order, neg_prot)

pheatmap(beta_mat[, pos_prot], display_numbers=text_mat[, pos_prot], 
  cluster_cols=FALSE, cluster_rows=FALSE,
  color=colorRampPalette(pal)(255), breaks=seq(-0.12, 0.12, length=256),
  legend_breaks=seq(-0.12, 0.10, by=0.02), border_color=NA,
  fontsize_row=8, fontsize_col=8, fontsize_number=8,
  cellwidth=8.6, cellheight=10,
  filename=sprintf("%s/sig_heatmap_pos.pdf", out_dir))

pheatmap(beta_mat[, neg_prot], display_numbers=text_mat[, neg_prot], 
  cluster_cols=FALSE, cluster_rows=FALSE,
  color=colorRampPalette(pal)(255), breaks=seq(-0.12, 0.12, length=256),
  legend_breaks=seq(-0.12, 0.10, by=0.02), border_color=NA,
  fontsize_row=8, fontsize_col=8, fontsize_number=8,
  cellwidth=8.6, cellheight=10,
  filename=sprintf("%s/sig_heatmap_neg.pdf", out_dir))


# -----------------------------------------------------------
# Validate any associations we can on olink panel
# -----------------------------------------------------------

olink_assocs <- foreach(grs = GRSs, .combine=rbind) %do% {
  dt <- fread(sprintf("analyses/univariate_associations/%s/olink_proteins/associations.tsv", grs))
  dt[, PRS := grs]
  dt
}
olink_info <- fread("analyses/processed_traits/olink_proteins/trait_info.tsv")
olink_assocs <- olink_assocs[olink_info, on = .(trait=variable), nomatch=0]
olink_assocs[, panel := gsub("_.*", "", trait)]
olink_assocs[grs_info, on = .(PRS=GRS_name), PRS := Display_name]
olink_assocs <- olink_assocs[panel != "neu"] # don't have agreements in place to use this data.

soma_assocs <- comb_assocs[Prot.FDR < 0.1, .(PRS, Gene, UniProt, 
                           Beta=Prot.Beta, L95=Prot.L95, 
                           U95=Prot.U95, Pvalue=Prot.Pvalue, FDR=Prot.FDR)]
soma_assocs[Gene == "PDE4D", c("Gene", "UniProt") := .("PDE4D|PDE4A", "Q08499|P27815")]
soma_assocs <- soma_assocs[, .(Gene = strsplit(Gene, "\\|")[[1]], UniProt=strsplit(UniProt, "\\|")[[1]]),
                           by=.(PRS, Beta, L95, U95, Pvalue, FDR)]

comp <- merge(soma_assocs, olink_assocs, by=c("PRS", "UniProt"))
comp <- comp[, .(PRS, UniProt, Gene,  
                 beta.soma=Beta, l95.soma=L95, u95.soma=U95, p.soma=Pvalue, fdr.soma=FDR,
                 beta.oli=beta, l95.oli=l95, u95.oli=u95, p.oli=pval, fdr.oli=fdr)]

# data.table for ribbon
plot_min <- comp[, min(c(l95.soma, l95.oli))]
plot_max <- comp[, max(c(u95.soma, u95.oli))]
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

g <- ggplot(comp) +
  aes(x = beta.soma, y = beta.oli,
      xmin = l95.soma, ymin = l95.oli,
      xmax = u95.soma, ymax = u95.oli,
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
  geom_text(inherit.aes=FALSE, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1,
            label = sprintf("Pearson r = %s", comp[,round(cor(beta.soma, beta.oli)*100)/100])) +
  scale_x_continuous(name = "Beta (95% CI) somalogic platform", expand=c(0,0)) +
  scale_y_continuous(name = "Beta (95% CI) olink platform", expand=c(0,0)) +
  scale_color_manual(name = "PRS", values=col_map, 
                     guide = guide_legend(direction="horizontal", nrow=3)) +
  scale_fill_manual(name = "PRS", values=col_map,
                    guide = guide_legend(direction="horizontal", nrow=3)) +
  coord_fixed(ratio=1) +
  theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g, width=7.2/2, height=7.2/2, file=sprintf("%s/soma_olink_compare.pdf", out_dir), useDingbats=FALSE)

# Output information for supp table
olink_prot <- fread("analyses/processed_traits/olink_proteins/traits.tsv", colClasses=c("IID"="character"))
olink_prot <- olink_prot[olink_info, on = .(variable), nomatch=0]
olink_prot[, panel := gsub("_.*", "", variable)]
olink_prot <- olink_prot[UniProt %in% comp$UniProt]
olink_prot <- olink_prot[, .(N=length(intersect(IID, levels$IID)), 
                             overlap=length(intersect(IID, prot$IID))),
                         by=.(UniProt, panel)]

olink_assocs <- comp[olink_prot, on = .(UniProt)]
olink_assocs <- olink_assocs[, .(PRS, UniProt, Gene, Panel=panel, N, Overlap=overlap,
																 Beta=beta.oli, L95=l95.oli, U95=u95.oli,
																 Pvalue=p.oli)]
olink_assocs <- olink_assocs[order(Pvalue)]
fwrite(olink_assocs, sep="\t", quote=FALSE, file=sprintf("%s/olink_assocs.tsv", out_dir))

# -----------------------------------------------------------
# Compare all olink assocaitions where we can
# -----------------------------------------------------------

olink_assocs <- foreach(grs = GRSs, .combine=rbind) %do% {
  dt <- fread(sprintf("analyses/univariate_associations/%s/olink_proteins/associations.tsv", grs))
  dt[, PRS := grs]
  dt
}
olink_info <- fread("analyses/processed_traits/olink_proteins/trait_info.tsv")
olink_assocs <- olink_assocs[olink_info, on = .(trait=variable), nomatch=0]
olink_assocs[, panel := gsub("_.*", "", trait)]
olink_assocs[grs_info, on = .(PRS=GRS_name), PRS := Display_name]

# Drop neurology panel
olink_assocs <- olink_assocs[panel != "neu"]

olink_assocs <- olink_assocs[, .(beta=mean(beta), l95=mean(l95), u95=mean(u95), pval=mean(pval)),
                             by=.(PRS, UniProt)]
olink_assocs[,fdr := p.adjust(pval, method="fdr"), by=PRS]
olink_assocs[olink_info, on = .(UniProt), Target := protein]


soma_assocs <- comb_assocs[, .(PRS, Target, UniProt, beta=Prot.Beta, 
                               l95=Prot.L95, u95=Prot.U95, pval=Prot.Pvalue, 
                               fdr=Prot.FDR)]
soma_assocs <- unique(soma_assocs)

comb <- merge(olink_assocs, soma_assocs, by = c("PRS", "UniProt"), suffixes=c(".olink", ".soma"))
comb[, group := "No association"]
comb[fdr.soma < 0.05 & fdr.olink >= 0.05, group := "FDR < 0.05 for SomaLogic"]
comb[fdr.soma >= 0.05 & fdr.olink < 0.05, group := "FDR < 0.05 for Olink"]
comb[fdr.soma < 0.05 & fdr.olink < 0.05, group := "FDR < 0.05 for both"]

g <- ggplot(comb, aes(y = beta.soma, x = beta.olink,
                      ymin = l95.soma, ymax = u95.soma,
                      xmin = l95.olink, xmax = u95.olink,
                      color=group, fill=group)) +
     geom_hline(yintercept=0, linetype=2) +
     geom_vline(xintercept=0, linetype=2) +
     geom_errorbarh(data=comb[group == "No association"], alpha=0.1) +
     geom_errorbar(data=comb[group == "No association"], alpha=0.1) +
     geom_point(data=comb[group == "No association"], pch=19, fill="black", alpha=0.3, size=0.3) +
     geom_errorbarh(data=comb[group != "No association"], alpha=0.5) +
     geom_errorbar(data=comb[group != "No association"], alpha=0.5) +
     geom_point(data=comb[group != "No association"], pch=21, color="black", size=0.8) +
     scale_y_continuous(name="Association with SomaLogic", breaks=c(-0.1, 0, 0.1)) +
     scale_x_continuous(name="Association with Olink assay", breaks=c(-0.1, 0, 0.1)) +
     scale_fill_manual(name="", values=c("No association"="black", "FDR < 0.05 for SomaLogic"="#e41a1c",
                                "FDR < 0.05 for Olink"="#ff7f00", "FDR < 0.05 for both"="#984ea3")) +
     scale_color_manual(name="", values=c("No association"="black", "FDR < 0.05 for SomaLogic"="#e41a1c",
                                 "FDR < 0.05 for Olink"="#ff7f00", "FDR < 0.05 for both"="#984ea3")) +
     facet_wrap(~ PRS, nrow=1) +
     ggtitle("Comparison of associations for 271 proteins measured by both Olink and SomaLogic") +
     theme_bw() + theme(
       axis.title=element_text(size=8, face=1), axis.text=element_text(size=8),
       strip.text=element_text(size=8), title=element_text(size=8, face=2),
       panel.grid=element_blank(), legend.position="bottom"
     )
ggsave(g, file=sprintf("%s/olink_soma_all_compare.pdf", out_dir), width=7.2, height=2.6, useDingbats=FALSE)

# -----------------------------------------------------------
# Sensitivity analysis to BMI, season, and circadian effects
# -----------------------------------------------------------

# Load associations
type <- c("BMI", "circadian", "season")
adj <- foreach(grs = GRSs, .combine=rbind) %do% {
  foreach(adj_type = type, .combine=rbind) %do% {
    dt <- fread(sprintf("analyses/sensitivity_associations/%s_adjusted/%s/somalogic_proteins/associations.tsv", adj_type, grs))
    dt[, c("grs", "model") := .(grs, adj_type)]
    dt
  }
}

# Merge information and collapse to protein level
adj[prot_info, on = .(trait=variable), Aptamer := SeqId]
adj[grs_info, on = .(grs = GRS_name), PRS := Display_name]
comp <- comb_assocs[Prot.FDR < 0.05]
adj <- adj[,.(PRS, Adjustment = model, Aptamer, Apt.Beta.Adj=beta, 
              Apt.L95.Adj=l95, Apt.U95.Adj=u95, Apt.Pvalue.Adj=pval)]
comp <- comp[adj, on = .(PRS, Aptamer), nomatch=0]
comp <- comp[, .(Prot.Beta.Adj=mean(Apt.Beta.Adj),
                 Prot.L95.Adj=mean(Apt.L95.Adj),
                 Prot.U95.Adj=mean(Apt.U95.Adj),
                 Prot.Pvalue.Adj=prot_pvalue(Apt.Pvalue.Adj, Apt.Beta.Adj)),
             by=.(PRS, Gene, UniProt, Prot.Beta, Prot.L95, Prot.U95, 
                  Prot.Pvalue, Prot.FDR, Adjustment)]

# data.table for ribbon
rdt <- comp[, .(x=c(min(c(Prot.L95, Prot.L95.Adj)), 0, 
                    max(c(Prot.U95, Prot.U95.Adj)))*1.05,
                ymax=c(min(c(Prot.L95, Prot.L95.Adj)), 0, 
                       max(c(Prot.U95, Prot.U95.Adj)))*1.05, 
                ymin=c(0,0,0)),
             by=PRS]
 
g <- ggplot(comp) +
	aes(x = Prot.Beta, y = Prot.Beta.Adj,
		  xmin = Prot.L95, ymin = Prot.L95.Adj,
		  xmax = Prot.U95, ymax = Prot.U95.Adj, 
      color = Adjustment, fill = Adjustment) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
	geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
	geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
	geom_errorbarh(height=0, alpha=0.7, size=0.5) +
	geom_errorbar(width=0, alpha=0.7, size=0.5) +
	geom_abline(intercept = 0, slope = 1, linetype=2, color="#737373") +
	geom_point(shape = 21, size=1, color="#00000000") +
	facet_grid(Adjustment ~ PRS, scale="free") +
	scale_x_continuous(name = "Beta (95% CI)", expand=c(0,0)) +
	scale_y_continuous(name = "Beta (95% CI)") +
  scale_color_manual(guide=FALSE, values=c("BMI"="#fdbb84", "season"="#addd8e", "circadian"="#41b6c4")) +
  scale_fill_manual(guide=FALSE, values=c("BMI"="#b30000", "season"="#238b45", "circadian"="#225ea8")) +
	theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank()
  )
ggsave(g, width=7.2, height=7.2, file=sprintf("%s/sensitivity.pdf", out_dir), useDingbats=FALSE)

# -----------------------------------------------------------
# Build BMI "mediation analysis" tables
# -----------------------------------------------------------

# Get information about the proteins ablated by bmi adjustment
bmi_prot <- comp[Prot.FDR < 0.1 & Prot.Pvalue.Adj >= 0.05 & Adjustment == "BMI", .(Gene, UniProt)]
bmi_prot <- prot_info[bmi_prot, on=.(Gene.Name=Gene, UniProt.Id.Current.at.Uniprot=UniProt), 
                      .(variable, SeqId, Gene, UniProt)]

# Get their aptamer levels
bmi_apts <- prot[bmi_prot, on = .(variable)]
bmi_apts <- dcast(bmi_apts, IID ~ variable, value.var="value")

# Build a data.table containing the protein levels, PRS, and relevant covariats
prot_covar <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv", colClasses=c("IID"="character"))
dat <- pheno[, .(IID, age=agePulse, sex=sexPulse, bmi)]
dat <- dat[prot_covar, on = .(IID), nomatch=0]
dat <- dat[levels[grs == "T2D_2018", .(IID, T2D_PRS=score)], on = .(IID), nomatch=0]
dat <- dat[pcs, on = .(IID), nomatch=0]
dat <- dat[bmi_apts, on = .(IID), nomatch=0]

# Function that gives the beta, l95, u95, and pvalue for linear regression
# 'vars' - model variables to extract
# '...' - arguments to pass to lm()
lm_coef <- function(vars, ...) { 
  l1 <- lm(...)
  cf <- coef(summary(l1))
  ci <- confint(l1)
  data.table(vars=vars, beta=cf[vars,1], l95=ci[vars,1], u95=ci[vars,2], pval=cf[vars,4])
} 

# First, refit the Protein ~ PRS + BMI model to also extract BMI coefficients
model <- paste("%s ~ T2D_PRS + %s", "age + factor(sex) + factor(batch)",
               "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10",
               sep=" + ")

# Run for each protein aptamer in turn
bmi_adj <- foreach(apt = bmi_prot$variable) %do% {
  lm_coef(c("T2D_PRS", "bmi"), sprintf(model, apt, "bmi"), data=dat)
}
names(bmi_adj) <- bmi_prot$variable
bmi_adj <- rbindlist(bmi_adj, idcol="aptamer")

# Now fit multivariable models BMI ~ PRS + protein
bmi_med <- foreach(apt = bmi_prot$variable) %do% {
    lm_coef(c("T2D_PRS", apt), sprintf(model, "bmi", apt), data=dat)
}
names(bmi_med) <- bmi_prot$variable
bmi_med <- rbindlist(bmi_med, idcol="aptamer")
bmi_med[vars != "T2D_PRS", vars := "aptamer"]

# Combine results into single table
bmi_dt <- rbind("BMI_adj"=bmi_adj, "BMI_med"=bmi_med, idcol="model")

# Add back protein information, then average across aptamers (multiple for WFIKKN2 only)
bmi_dt <- bmi_dt[bmi_prot, on = .(aptamer=variable)]
bmi_dt <- bmi_dt[, .(beta=mean(beta), l95=mean(l95), u95=mean(u95), pval=mean(pval)), 
                 by=.(model, Gene, UniProt, vars)]

# Cast out to wide format
bmi_dt <- dcast(bmi_dt, Gene + UniProt ~ model + vars, value.var=c("beta", "l95", "u95", "pval"))

# Add in univariate unadjusted information for reference:
bmi_unadj <- comp[Prot.FDR < 0.1 & Prot.Pvalue.Adj >= 0.05 & Adjustment == "BMI", 
                  .(Gene, UniProt, Group = ifelse(Prot.FDR < 0.05, "FDR < 0.05", "FDR < 0.1"),
                    beta_unadj=Prot.Beta, l95_unadj=Prot.L95, u95_unadj=Prot.U95, 
                    pval_unadj=Prot.Pvalue, fdr_unadj=Prot.FDR)]
bmi_dt <- bmi_unadj[bmi_dt, on = .(Gene, UniProt)]

# Order rows and columns and write out
bmi_dt <- bmi_dt[order(-pval_BMI_adj_T2D_PRS)][order(Group)]
bmi_dt <- bmi_dt[, .(Protein=Gene, UniProt, Group, 
                     beta_unadj, l95_unadj, u95_unadj, pval_unadj, fdr_unadj,
                     beta_BMI_adj_T2D_PRS, l95_BMI_adj_T2D_PRS, u95_BMI_adj_T2D_PRS, pval_BMI_adj_T2D_PRS,
                     beta_BMI_adj_bmi, l95_BMI_adj_bmi, u95_BMI_adj_bmi, pval_BMI_adj_bmi,
                     beta_BMI_med_T2D_PRS, l95_BMI_med_T2D_PRS, u95_BMI_med_T2D_PRS, pval_BMI_med_T2D_PRS,
                     beta_BMI_med_aptamer, l95_BMI_med_aptamer, u95_BMI_med_aptamer, pval_BMI_med_aptamer)]
fwrite(bmi_dt, sep="\t", quote=FALSE, file=sprintf("%s/bmi_moderated.tsv", out_dir))

# -----------------------------------------------------------
# Multivariable models of proteins associated with > 1 PRS
# -----------------------------------------------------------

mult <- comb_assocs[Prot.FDR < 0.1, .N, by=.(Gene, PRS)][,.N,by=Gene][N > 1]
mult <- comb_assocs[Prot.FDR < 0.1][mult, on = .(Gene)]
prot_pairs <- unique(mult[,.(Gene, Aptamer, PRS)])
prot_pairs[prot_info, on = .(Aptamer=SeqId), variable := variable]
prot_pairs[grs_info, on = .(PRS=Display_name), grs := GRS_name]

covar <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv", colClasses=c("IID"="character"))
covar <- covar[pcs, on = .(IID), nomatch=0]

# Fit multivariable models for each aptamer 
multi_prs <- foreach(aptvar = unique(prot_pairs$variable), .combine=rbind) %do% {
  prss <- prot_pairs[variable == aptvar, grs]

  test_prs_levels <- dcast(levels[grs %in% prss], IID ~ grs, value.var="score")
  test_apt_levels <- dcast(prot[variable == aptvar], IID ~ variable, value.var="value")
  dat <- test_apt_levels[test_prs_levels, on =.(IID), nomatch=0][covar, on = .(IID), nomatch=0]
  model <- sprintf("%s ~ %s + factor(batch)", aptvar, paste(prss, collapse=" + "))
  model <- paste(model, "+", paste(paste0("PC_", 1:10), collapse=" + "))

  l1 <- lm(as.formula(model), data=dat)
  l1_ci <- confint(l1)
  l1_coef <- coef(summary(l1))
  data.table(PRS = prot_pairs[variable == aptvar, PRS],
             Gene = prot_pairs[variable == aptvar, Gene],
             Aptamer = prot_pairs[variable == aptvar, Aptamer],
             beta=l1_coef[prss, 1], l95=l1_ci[prss, 1], 
             u95=l1_ci[prss,2], pval=l1_coef[prss, 4])
}

# Collapse to protein level (relevant for SHBG only)
multi_prs <- multi_prs[, .(Prot.Beta = mean(beta), Prot.L95 = mean(l95), 
                           Prot.U95 = mean(u95), Prot.Pvalue = prot_pvalue(pval, beta)),
                        by=.(PRS, Gene)]

comp <- merge(multi_prs, unique(comb_assocs, by=c("PRS", "Gene")), 
              by=c("Gene", "PRS"), suffixes=c(".Multi", ""))

# data.table for ribbon
rdt <- comp[, .(x=c(min(c(Prot.L95, Prot.L95.Multi)), 0, 
                    max(c(Prot.U95, Prot.U95.Multi)))*1.05,
                ymax=c(min(c(Prot.L95, Prot.L95.Multi)), 0, 
                       max(c(Prot.U95, Prot.U95.Multi)))*1.05, 
                ymin=c(0,0,0)),
             by=PRS]

g <- ggplot(comp) +
  aes(x = Prot.Beta, y = Prot.Beta.Multi,
      xmin = Prot.L95, ymin = Prot.L95.Multi,
      xmax = Prot.U95, ymax = Prot.U95.Multi,
      color=PRS, fill=PRS) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
	geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
	geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#737373") +
  geom_errorbarh(height=0, alpha=1, size=0.5) +
  geom_errorbar(width=0, alpha=1, size=0.5) +
  geom_point(shape = 21, size=1.3, color="#00000000") +
  scale_x_continuous(name = "Beta (95% CI) when testing each PRS separately", expand=c(0,0)) +
  scale_y_continuous(name = "Beta (95% CI) for PRSs in multivariable model", expand=c(0,0)) +
  scale_color_manual(name = "PRS", values=c("Coronary Artery Disease"="#b2abd2", "Type 2 Diabetes"="#fdb863")) +
  scale_fill_manual(name = "PRS", values=c("Coronary Artery Disease"="#5e3c99", "Type 2 Diabetes"="#e66101")) +
  facet_wrap(~ Gene, nrow=1) +
  coord_fixed(ratio=1) +
  theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank(), 
    legend.position="bottom"
  )
ggsave(g, width=7.2, height=3.15, file=sprintf("%s/multi_prs.pdf", out_dir), useDingbats=FALSE)


# -----------------------------------------------------------
# Sensitivity analysis to pQTL adjustment
# -----------------------------------------------------------

pqtl_adj <- foreach(grs = GRSs) %do% {
  fread(sprintf("analyses/sensitivity_associations/qtl_prob_dosage_adjusted/%s/somalogic_proteins/qtl_grs_associations.tsv", grs))
}
names(pqtl_adj) <- GRSs
pqtl_adj <- rbindlist(pqtl_adj, idcol="grs")
pqtl_adj[grs_info, on = .(grs=GRS_name), PRS := Display_name]
pqtl_adj[prot_info, on = .(trait = variable), Aptamer := SeqId]

# Filter
pqtl_adj <- pqtl_adj[comb_assocs[Prot.FDR < 0.1, .(PRS, Aptamer, Gene, UniProt)], on = .(PRS, Aptamer), nomatch=0]

# Flag proteins that actually have pQTLs
has_pQTLs <- unique(pqtl_adj[variable != "grs", .(PRS, Gene, Aptamer)])

# Summarise at the protein level
pqtl_prot_adj <- pqtl_adj[variable == "grs", 
                          .(Prot.Beta.Adj=mean(beta), Prot.L95.Adj=mean(l95),
                            Prot.U95.Adj=mean(u95), Prot.Pvalue.Adj=prot_pvalue(pval, beta)),
                          by=.(PRS, Gene, UniProt)]

# Combine tables
pqtl_adj <- pqtl_adj[,.(PRS, Gene, UniProt, Aptamer, variable, Apt.Beta=beta, 
                        Apt.L95=l95, Apt.U95=u95, Apt.Pvalue=pval)]
pqtl_adj <- pqtl_prot_adj[pqtl_adj, on = .(PRS, Gene, UniProt)]


# Compare
comp <- merge(pqtl_prot_adj, unique(comb_assocs[Prot.FDR < 0.1], by=c("PRS", "Gene", "UniProt")),
              by=c("PRS", "Gene", "UniProt"), all.y=TRUE)

# Annotate
comp[, type := "has_pQTLs"]
comp[!has_pQTLs, on = .(PRS, Gene), type := "no_pQTLs"]

# Order points
comp[, type := factor(type, levels=c("no_pQTLs", "has_pQTLs"))]
comp <- comp[order(type)]

# data.table for ribbon
rdt <- comp[, .(x=c(min(c(Prot.L95, Prot.L95.Adj)), 0, 
                    max(c(Prot.U95, Prot.U95.Adj)))*1.05,
                ymax=c(min(c(Prot.L95, Prot.L95.Adj)), 0, 
                       max(c(Prot.U95, Prot.U95.Adj)))*1.05, 
                ymin=c(0,0,0)),
             by=PRS]

g <- ggplot(comp) +
	aes(x = Prot.Beta, y = Prot.Beta.Adj,
		  xmin = Prot.L95, ymin = Prot.L95.Adj,
		  xmax = Prot.U95, ymax = Prot.U95.Adj, 
      color = type, fill = type) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
	geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
	geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
	geom_errorbarh(height=0, alpha=0.7, size=0.5) +
	geom_errorbar(width=0, alpha=0.7, size=0.5) +
	geom_abline(intercept = 0, slope = 1, linetype=2, color="#bdbdbd") +
	geom_point(shape = 21, size=1.2, color="#00000000") +
	facet_wrap(~ PRS, nrow=1, scales="free_x") +
	scale_x_continuous(name = "Beta (95% CI) for the PRS", expand=c(0,0)) +
	scale_y_continuous(name = "Beta (95% CI) adjusting for conditionally independent pQTLs") +
  scale_fill_manual(guide=FALSE, 
    values=c("has_pQTLs"="#1f78b4", "no_pQTLs"="#525252")) +
  scale_color_manual(guide=FALSE, 
    values=c("has_pQTLs"="#a6cee3", "no_pQTLs"="#969696")) +
	theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank()
  )
ggsave(g, width=7.2, height=2.8, file=sprintf("%s/pqtl_adj.pdf", out_dir), useDingbats=FALSE)

# For paper: how many FDR < 0.05 proteins have at least one pQTL?
pqtl_adj[has_pQTLs, on = .(PRS, Gene)][
         comb_assocs[Prot.FDR < 0.05, .(PRS, Gene)], on = .(PRS, Gene), nomatch=0][
         variable == "grs", length(unique(Gene))]
# For paper: how many FDR < 0.1 proteins have at least one pQTL?
pqtl_adj[has_pQTLs, on = .(PRS, Gene)][
         comb_assocs[Prot.FDR < 0.1, .(PRS, Gene)], on = .(PRS, Gene), nomatch=0][
         variable == "grs", length(unique(Gene))]

# Tabulate contributions - also output information about each relevant pQTL
cond_pQTLs <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=5, startRow=3)
cond_pQTLs <- as.data.table(cond_pQTLs)

head_row <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, rows=5:6, fillMergedCells=TRUE)
head_row <- gsub(".NA$", "", paste(colnames(head_row), as.vector(head_row[1,]), sep="."))
pQTLs <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, startRow=6)
colnames(pQTLs) <- head_row
pQTLs <- as.data.table(pQTLs)
cis_pQTLs <- pQTLs[`cis/.trans` == "cis"]

# Summarise information into a single string
pqtl_info <- pQTLs[cond_pQTLs, on=.(SOMAmer.ID, `Sentinel.variant*`=Sentinel.variant)]
pqtl_info <- pqtl_info[, .(SOMAmer.ID, pQTL=Conditional.variant, sentinel=`Sentinel.variant*`, 
                           Chr, Pos, sentinel_ld=`LD.with.sentinel.variant.(r2)`,
                           cis_or_trans=`cis/.trans`)]
pqtl_info[, ld_text := paste0("r2=", round(sentinel_ld*100)/100)]
pqtl_info[pQTL == sentinel, msg := paste0(cis_or_trans, "-pQTL. Sentinel variant on chromosome ", Chr, 
                                          " at ", format(Pos, big.mark=",", trim=TRUE), ".")]
pqtl_info[pQTL != sentinel, msg := paste0(cis_or_trans, "-pQTL conditionally independent of ", sentinel,
                                          " (", ld_text, ").")]
pqtl_info <- pqtl_info[, .(SOMAMER_ID=SOMAmer.ID, pQTL, pqtl_info=msg)]
pqtl_info[prot_info, on = .(SOMAMER_ID), Aptamer := SeqId]

# Create output table - first pull in all required information
pqtl_adj <- pqtl_adj[has_pQTLs, on = .(PRS, Gene, Aptamer)]
pqtl_adj[comb_assocs[Prot.FDR < 0.05], on = .(PRS, Gene), type := "FDR < 0.05"]
pqtl_adj[is.na(type), type := "FDR < 0.1"]
pqtl_adj[pqtl_info, on = .(Aptamer, variable=pQTL), pQTL_info := i.pqtl_info]

# Filter to relevant columns
pqtl_adj <- pqtl_adj[, .(PRS, Protein=Gene, UniProt, Group=type, 
                       Prot.Beta.Adj, Prot.L95.Adj, Prot.U95.Adj, Prot.Pvalue.Adj,
                       Aptamer, variable, Apt.Beta, Apt.L95, Apt.U95, Apt.Pvalue,
                       pQTL_info)]
pqtl_adj[variable == "grs", variable := "PRS"]

# Order rows:
pqtl_adj <- pqtl_adj[order(Apt.Pvalue)]
pqtl_adj <- pqtl_adj[order(-grepl("Sentinel", pQTL_info))]
pqtl_adj <- pqtl_adj[order(-grepl("cis", pQTL_info))]
pqtl_adj <- pqtl_adj[order(variable != "PRS")]
pqtl_adj <- pqtl_adj[order(Aptamer)]
pqtl_adj <- pqtl_adj[order(Protein)]
pqtl_adj <- pqtl_adj[order(Prot.Pvalue.Adj)]
pqtl_adj <- pqtl_adj[order(Group)]
pqtl_adj <- pqtl_adj[order(PRS)]

# Determine effect allele
dosage <- fread("analyses/processed_traits/somalogic_proteins/qtl_geno_prob.tsv")
dosage <- melt(dosage, id.vars=c("IID"), variable.name="rsid", value.name="dosage")
alleles <- fread("analyses/processed_traits/somalogic_proteins/qtl_alleles.tsv")
alleles <- melt(alleles, id.vars=c("IID"), variable.name="rsid", value.name="alleles")
dosage_alleles <- merge(dosage, alleles, by=c("IID", "rsid"))
dosage_alleles <- dosage_alleles[rsid %in% pqtl_adj$variable]
rm(dosage, alleles); gc() # free up some memory

dosage_alleles[, A1 := sapply(strsplit(alleles, "/"), `[`, 1)]
dosage_alleles[, A2 := sapply(strsplit(alleles, "/"), `[`, 2)]
dosage_alleles <- dosage_alleles[A1 == A2]
dosage_alleles <- unique(dosage_alleles, by=c("rsid", "A1", "A2"))
dosage_alleles <- dosage_alleles[alleles != "0/0"]
dosage_alleles <- dosage_alleles[,.SD[which.max(dosage)],by=rsid]
dosage_alleles <- dosage_alleles[, .(rsid, EA=A1)]

pqtl_adj[dosage_alleles, on = .(variable = rsid), variable := sprintf("%s (%s)", variable, EA)]

fwrite(pqtl_adj, quote=FALSE, sep="\t", file=sprintf("%s/pqtl_adj.tsv", out_dir))

# -----------------------------------------------------------
# Sensitivity analysis to pQTL removal
# -----------------------------------------------------------

# Load associations after pQTL removal
sig_pairs <- comb_assocs[Prot.FDR < 0.1, .(PRS, Aptamer)]
sig_pairs[prot_info, on=.(Aptamer=SeqId), variable := variable]
sig_pairs[grs_info, on = .(PRS=Display_name), grs := GRS_name]
pqtl_removed <- foreach(ii = seq_len(nrow(sig_pairs)), .combine=rbind, .errorhandling="remove") %do% {
  sig_pairs[ii, fread(file=sprintf("analyses/grs_pqtl_removed/%s/%s/associations.tsv", grs, variable))]
}
pqtl_removed[grs_info, on = .(grs=GRS_name), PRS := Display_name]
pqtl_removed[prot_info, on = .(trait=variable), Aptamer := SeqId]
pqtl_removed[comb_assocs, on = .(PRS, Aptamer), c("Gene", "UniProt") := .(Gene, UniProt)]

# collapse across proteins
pqtl_removed <- pqtl_removed[, .(Prot.Beta.Adj = mean(beta), Prot.L95.Adj = mean(l95),
                                 Prot.U95.Adj = mean(u95), Prot.Pvalue.Adj = mean(pval)),
                              by=.(PRS, Gene, UniProt)]

# Compare associations
comp <- merge(pqtl_removed, comb_assocs, by=c("PRS", "Gene", "UniProt"))

# data.table for ribbon
rdt <- comp[, .(x=c(min(c(Prot.L95, Prot.L95.Adj)), 0, 
                    max(c(Prot.U95, Prot.U95.Adj)))*1.05,
                ymax=c(min(c(Prot.L95, Prot.L95.Adj)), 0, 
                       max(c(Prot.U95, Prot.U95.Adj)))*1.05, 
                ymin=c(0,0,0)),
             by=PRS]

g <- ggplot(comp) +
	aes(x = Prot.Beta, y = Prot.Beta.Adj,
		  xmin = Prot.L95, ymin = Prot.L95.Adj,
		  xmax = Prot.U95, ymax = Prot.U95.Adj) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
  geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.7, size=0.5, color="#addd8e") +
  geom_errorbar(width=0, alpha=0.7, size=0.5, color="#addd8e") +
	geom_abline(intercept = 0, slope = 1, linetype=2, color="#bdbdbd") +
	geom_point(shape = 21, size=1.2, color="#00000000", fill="#238b45") +
	facet_wrap(~ PRS, nrow=1, scales = "free_x") +
	scale_x_continuous(name = "Beta (95% CI) for the PRS", expand=c(0,0)) +
	scale_y_continuous(name = "Beta (95% CI) removing pQTLs") +
	theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank()
  )
ggsave(g, width=7.2, height=2.8, file=sprintf("%s/pqtl_removed.pdf", out_dir), useDingbats=FALSE)

# --------------------------------------------------------------
# Collate MR and coloc results
# --------------------------------------------------------------

# Load results
instruments <- foreach(grs = GRSs, .combine=rbind, .errorhandling="remove") %do% {
  fread(sprintf("analyses/mendelian_randomisation/%s/instruments.tsv", grs))
}

mr <- foreach(grs = GRSs, .combine=rbind, .errorhandling="remove") %do% {
  fread(sprintf("analyses/mendelian_randomisation/%s/mr_results.tsv", grs))
} 

# Add PRS display names
mr[grs_info, on = .(PRS=GRS_name), PRS := Display_name]
instruments[grs_info, on = .(PRS=GRS_name), PRS := Display_name]

# Require proteins to have at least three pQTLs, at least one cis
prot_ivs <- unique(instruments, by=c("PRS", "GWAS", "Gene", "chr", "pos", "EA", "OA"))
iv_stats <- prot_ivs[!is.na(P.gwas), .(N=.N, cis=sum(type == "cis")), by=.(PRS, GWAS, Gene)]
pass <- iv_stats[N >= 3 & cis > 0, .(PRS, GWAS, Gene)]

# Drop ApoE - any causal effect likely reflects its 
# cross-reactivity to the ApoE2 and ApoE4 isoform aptamers,
# which we know will have a causal effect due to their 
# effects on LDL levels
pass <- pass[Gene != "APOE"]

# Filter the MR and IV tables
instruments <- instruments[pass, on = .(PRS, GWAS, Gene)]
prot_ivs <- prot_ivs[pass, on = .(PRS, GWAS, Gene)]
mr <- mr[pass, on = .(PRS, GWAS, Gene)]

# Restrict to representative subset of methods
mr <- rbind(
  mr[method %in% c("IVW", "Simple median", "Weighted median", "Weighted mode (simple SE)")],
  mr[c(which(method == "MR-Egger"), which(method == "MR-Egger") + 1)]
)
mr[method == "Weighted mode (simple SE)", method := "Weighted mode"]
mr <- mr[order(Gene)][order(PRS)]

# Drop methods that could not be run
mr <- mr[!is.na(mr_estimate)]

# Add instrument stats to mr table
mr <- mr[iv_stats, on = .(PRS, GWAS, Gene), nomatch=0]

# Transform log odds to odds ratios
mr[, c("mr_estimate", "mr_L95", "mr_U95") := .(exp(mr_estimate), exp(mr_L95), exp(mr_U95))]

# Get median across estimators:
mr_summary <- mr[method != "(intercept)", 
                 .(med_estimate=median(mr_estimate),
                   med_L95=median(mr_L95),
                   med_U95=median(mr_U95),
                   med_pval=median(mr_pval)),
                by = .(PRS, GWAS, Target, UniProt, Gene, N, cis)]
mr_summary[mr[method == "(intercept)"], on = .(PRS, GWAS, Target, UniProt, Gene, N, cis), horiz_pleio := mr_pval]

# Add to MR table and output
mr <- merge(mr_summary, mr, by = c("PRS", "GWAS", "Target", "UniProt", "Gene", "N", "cis"))
mr <- mr[order(med_pval)][order(horiz_pleio < 0.05)]
fwrite(mr, file=sprintf("%s/mr_results.tsv", out_dir), sep="\t", quote=FALSE)

# The plotting code is set up to work on log odds (centered at 0) not 
# odds ratios (centered at 1) due to the way the ribbon and polygon code
# is written, so we will convert back to log odds to set up the data.frames
# then back to odds ratios for the plots
mr[, c("mr_estimate", "mr_L95", "mr_U95") := .(log(mr_estimate), log(mr_L95), log(mr_U95))]
mr_summary[, c("med_estimate", "med_L95", "med_U95") := .(log(med_estimate), log(med_L95), log(med_U95))]

# Plot dose response curves for significant causal estimates:
sig_mr <- mr_summary[med_pval < 0.05 & horiz_pleio >= 0.05]
sig_mr_ivs <- instruments[sig_mr, on = .(PRS, GWAS, Target, UniProt, Gene)]
sig_mr_ivs <- sig_mr_ivs[, .(PRS, GWAS, Target, UniProt, Gene, type,
                             Prot.Effect.pQTL, Prot.SE.pQTL, Prot.P.pQTL, 
                             effect.gwas, se.gwas, P.gwas)]
sig_mr_ivs <- unique(sig_mr_ivs)

# ribbon for causal estimates
rdt <- sig_mr[, .(
         x = c(-1, 0, 1),
         ymin = c(-med_L95, 0, med_L95),
         ymax = c(-med_U95, 0, med_U95)),
        by = .(PRS, GWAS, Target, UniProt, Gene)]
plotlim <- sig_mr_ivs[, .(x=c(-1, 0, 1), xmult = c(
                          abs(min(Prot.Effect.pQTL - Prot.SE.pQTL)), 0,
                          max(Prot.Effect.pQTL + Prot.SE.pQTL))),
                      by=.(PRS, GWAS, Target, UniProt, Gene)]
expand <- plotlim[, .(expand=(max(xmult) - min(xmult))*0.05), 
                  by=.(PRS, GWAS, Target, UniProt, Gene)]
plotlim[expand, on = .(PRS, GWAS, Target, UniProt, Gene), xmult := xmult + expand]
rdt[plotlim, on = .(PRS, GWAS, Target, UniProt, Gene, x), 
    c("x", "ymin", "ymax") := .(x * xmult, ymin * xmult, ymax * xmult)]

g <- ggplot(sig_mr_ivs) +
  aes(x = Prot.Effect.pQTL, y = exp(effect.gwas),
      xmin = Prot.Effect.pQTL - Prot.SE.pQTL, 
      ymin = exp(effect.gwas - se.gwas),
      xmax = Prot.Effect.pQTL + Prot.SE.pQTL,
      ymax = exp(effect.gwas + se.gwas),
      fill = type) +
  geom_ribbon(data=rdt, aes(x=x, ymin=exp(ymin), ymax=exp(ymax)), inherit.aes=FALSE, fill="#ffeda0") +
  geom_abline(data=sig_mr, aes(intercept=1, slope=-(1-exp(med_estimate))), linetype=2, color="#f03b20") +
  geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, size=0.5, color="black", alpha=0.7) +
  geom_errorbar(width=0, size=0.5, color="black", alpha=0.7) +
	geom_point(shape = 21, size=1, color="black") +
	facet_wrap(~ PRS + Gene, nrow=1, scales="free") +
	scale_x_continuous(name = "SD effect on protein", expand=c(0,0)) +
	scale_y_continuous(name = "Odds Ratio") +
  scale_fill_manual(values=c("trans"="#ffcc00", "cis"="black")) +
	theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g, width=7.2, height=3.4, file=sprintf("%s/sig_mr.pdf", out_dir), useDingbats=FALSE)

# Generate dose response curves for all tested protein to disease pairs

# Get instruments for each protein:
prot_ivs <- unique(instruments[!is.na(effect.gwas), 
                                 .(PRS, GWAS, Target, UniProt, Gene, type,
                                   Prot.Effect.pQTL, Prot.SE.pQTL, Prot.P.pQTL,
                                   effect.gwas, se.gwas, P.gwas)])

# Determine the plot limits for each panel so we can derive the polygons
plotlim <- prot_ivs[, .(xmin = min(Prot.Effect.pQTL - Prot.SE.pQTL),
                        xmax = max(Prot.Effect.pQTL + Prot.SE.pQTL),
                        ymin = min(effect.gwas - se.gwas),
                        ymax = max(effect.gwas + se.gwas)),
                    by=.(PRS, GWAS, Target, UniProt, Gene)]

# Make sure to include intercepts
plotlim[xmax < 0, xmax := 0]
plotlim[xmin > 0, xmin := 0]
plotlim[ymax < 0, ymax := 0]
plotlim[ymin > 0, ymin := 0]

# Expand by 5% as per ggplot defaults
expand <- plotlim[, .(xexpand = abs((xmax - xmin))*0.05,
                      yexpand = abs((ymax - ymin))*0.05),
                  by=.(PRS, GWAS, Target, UniProt, Gene)]
plotlim[expand, on = .(PRS, GWAS, Target, UniProt, Gene),
        c("xmin", "xmax", "ymin", "ymax") := 
        .(xmin - xexpand, xmax + xexpand, ymin - yexpand, ymax + yexpand)]

# In a few cases the ylimit is not expanded because ggplot has no knowledge of
# these limits above so constructs limits based on the points, lines, and polygons
# derived below. To make sure each facet follows the above, we construct a data.table
# of points at each of the four corners to plot invisibly.
facetlim <- plotlim[, .(x=c(xmin, xmax, xmax, xmin),
                        y=c(ymin, ymin, ymax, ymax)),
                    by = .(PRS, GWAS, Target, UniProt, Gene)]

# Build table of information for dose response/causal estimate lines
drc <- mr[method != "(intercept)", .(PRS, GWAS, Target, UniProt, Gene, method, intercept=0, mr_estimate)]
eg_intercept <- mr[method == "(intercept)", .(PRS, GWAS, Target, UniProt, Gene, 
                                              method="MR-Egger", intercept=mr_estimate)]
drc[eg_intercept, on = .(PRS, GWAS, Target, UniProt, Gene, method), intercept := i.intercept]
med <- mr_summary[, .(PRS, GWAS, Target, UniProt, Gene, method="Median Estimate", 
                      intercept=0, mr_estimate=med_estimate)]
drc <- rbind(drc, med)

# Now build polygons for 95% CIs. This table is a series of x and y-coordinates for drawing
# the polygon for each method and protein to disease pair. For each method, we need to work
# out where the polygon intersects the edges of the plot so we know where to define the 
# polygon edges. See src/07_job_scripts/07_helpers/mr_functions.R for function code:
ci95 <- foreach(testIdx = plotlim[,.I], .combine=rbind) %do% {
  foreach(meth = unique(drc$method), .combine=rbind) %do% {
    test <- plotlim[testIdx, .(PRS, GWAS, Target, UniProt, Gene)]
    mrtest <- mr[test, on = .(PRS, GWAS, Target, UniProt, Gene)]
    if (meth == "MR-Egger") {
      intercept.l95 <- mrtest[method == "(intercept)", mr_L95]
      intercept.u95 <- mrtest[method == "(intercept)", mr_U95]
    } else {
      intercept.l95 <- 0
      intercept.u95 <- 0
    }

    if (meth == "Median Estimate") {
      slope.l95 <- mrtest[, unique(med_L95)]
      slope.u95 <- mrtest[, unique(med_U95)]
    } else {
      slope.l95 <- mrtest[method == meth, mr_L95]
      slope.u95 <- mrtest[method == meth, mr_U95]
    }

    xlim <- plotlim[testIdx, c(xmin, xmax)]
    ylim <- plotlim[testIdx, c(ymin, ymax)]

    if (intercept.l95 < ylim[1]) intercept.l95 <- ylim[1]
    if (intercept.u95 > ylim[2]) intercept.u95 <- ylim[2]

    poly <- line_ci95_poly(intercept.l95, intercept.u95, slope.l95, slope.u95, xlim, ylim)
    cbind(test, data.table(method=meth), poly)
  }
}

# Create short label for plot
short_prs <- data.table(Long = c("Chronic Kidney Disease", "Coronary Artery Disease", "Type 2 Diabetes"),
                        Short = c("CKD", "CAD", "T2D"))
prot_ivs[short_prs, on = .(PRS = Long), Label := sprintf("%s on %s", Gene, Short)]
ci95[short_prs, on = .(PRS = Long), Label := sprintf("%s on %s", Gene, Short)]
drc[short_prs, on = .(PRS = Long), Label := sprintf("%s on %s", Gene, Short)]

# Get order
mro <- unique(mr[,.(PRS, GWAS, Target, UniProt, Gene)])
mro[short_prs, on = .(PRS = Long), Label := sprintf("%s on %s", Gene, Short)]
facetlim[short_prs, on = .(PRS = Long), Label := sprintf("%s on %s", Gene, Short)]
prot_ivs[, Label := factor(Label, levels=mro$Label)]
ci95[, Label := factor(Label, levels=mro$Label)]
drc[, Label := factor(Label, levels=mro$Label)]
facetlim[, Label := factor(Label, levels=mro$Label)]

g <- ggplot(prot_ivs) +
  aes(x = Prot.Effect.pQTL, y = exp(effect.gwas),
      xmin = Prot.Effect.pQTL - Prot.SE.pQTL, 
      ymin = exp(effect.gwas - se.gwas),
      xmax = Prot.Effect.pQTL + Prot.SE.pQTL,
      ymax = exp(effect.gwas + se.gwas)) +
  geom_point(data=facetlim, inherit.aes=FALSE, aes(x=x, y=exp(y)), size=0) + # invisible points for determining facet limits.
  geom_polygon(data=ci95, aes(x=x, y=exp(y), fill=method, color=method), inherit.aes=FALSE, size=0.05) +
  scale_fill_manual(values=c("IVW"="#f8766d1a", "Simple median"="#a3a5001a", 
                             "Weighted median"="#7cae001a", "Weighted mode"="#00bfc41a",
                             "MR-Egger"="#c77cff1a", "Median Estimate"="#ffeda07d")) +
  scale_color_manual(values=c("IVW"="#f8766d66", "Simple median"="#a3a50066", 
                             "Weighted median"="#7cae0066", "Weighted mode"="#00bfc433",
                             "MR-Egger"="#c77cff66", "Median Estimate"="#f0d320c4")) +
  new_scale_fill() + new_scale_color() +
  geom_abline(data=drc, aes(intercept=exp(intercept), slope=-(1-exp(mr_estimate)), color=method), linetype=2) +
  scale_color_manual(values=c("IVW"="#f8766dff", "Simple median"="#a3a500ff", 
                             "Weighted median"="#7cae00ff", "Weighted mode"="#00bfc466",
                             "MR-Egger"="#c77cffff", "Median Estimate"="#f03b20ff")) +
  new_scale_color() +
  geom_hline(yintercept=1, color="#636363", linetype="dotted") +
  geom_vline(xintercept=0, color="#636363", linetype="dotted") +
  geom_errorbarh(height=0, size=0.5, color="black", alpha=0.7) +
  geom_errorbar(width=0, size=0.5, color="black", alpha=0.7) +
	geom_point(shape = 21, size=1, color="black", aes(fill=type)) +
	facet_wrap(~ Label, ncol=4, scales="free") +
	scale_x_continuous(name = "SD effect on protein", expand=c(0,0)) +
	scale_y_continuous(name = "Odds Ratio", expand=c(0,0)) +
  scale_fill_manual(values=c("trans"="#ffcc00", "cis"="black")) +
	theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g, width=7.2, height=8.3, file=sprintf("%s/all_mr.pdf", out_dir), useDingbats=FALSE)

# Format and output all instruments
instruments[type2 == "sub-threshold", type := "cis*"]
instruments[prot_info, on = .(SOMAMER_ID), Aptamer := SeqId]

instruments <- instruments[, .(PRS, Gene, UniProt, 
  rsID=RSID, chr, pos, EA, OA, "cis/trans"=type,
  sentinel, match, proxy_for, r2=r2_with_proxy, proxies, matched,
  EAF.pqtl, Prot.Effect.pQTL, Prot.SE.pQTL, Prot.P.pQTL,
  EAF.gwas, effect.gwas, se.gwas, P.gwas,
  Aptamer, aptamer_pqtl, Apt.Effect.pQTL, Apt.SE.pQTL, Apt.P.pQTL)]
mo <- unique(mr[,.(PRS, UniProt, Gene)])
instruments <- instruments[mo, on = .(PRS, UniProt, Gene)]

fwrite(instruments, file=sprintf("%s/mr_instruments.tsv", out_dir), sep="\t", quote=FALSE)

