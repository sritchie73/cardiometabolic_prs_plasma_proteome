library(data.table)
library(foreach)
library(ggplot2)
library(ggrastr)
library(cowplot)
library(scales)

# Test sensitivity of associations to PGS construction method

out_dir <- "analyses/pub/cardiometabolic_proteins/review1"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# load analysed PRSs
view <- fread("views/cardiometabolic.txt")

prs <- foreach(prs_id = view$GRS_name, .combine=rbind) %do% {
  this_prs <- fread(sprintf("analyses/GRS_profiles/%s/profile.sscore.gz", prs_id))
  this_prs[, PRS := prs_id]
  this_prs[, disease := view[GRS_name == prs_id, Display_name]]
  this_prs
}
prs <- prs[, .(IID, PRS, disease, score_sum)]

# Load PCs
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Filter to analysis samples
soma_prot <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")
soma_samples <- unique(soma_prot$IID)
rm(soma_prot); gc()

prs <- prs[IID %in% soma_samples]
pcs <- pcs[IID %in% soma_samples]

# Load in PRSs we want to compare to:
comp <- rbind(idcol="PRS",
  "Khera2018_T2D_PGS000014"=fread("analyses/GRS_profiles/Khera2018_T2D_PGS000014/profile.sscore.gz"),
  "Khera2018_CAD_PGS000013"=fread("analyses/GRS_profiles/Khera2018_CAD_PGS000013/profile.sscore.gz"),
  "Khera2018_afib_PGS000016"=fread("analyses/GRS_profiles/Khera2018_afib_PGS000016/profile.sscore.gz"),
  "Mahajan2018_T2D_PGS000036"=fread("analyses/GRS_profiles/Mahajan2018_T2D_PGS000036/profile.sscore.gz"),
  "Weng2017_afib_PGS000035"=fread("analyses/GRS_profiles/Weng2017_afib_PGS000035/profile.sscore.gz"),
  "Wuttke2019_eGFR_PGSuncurated"=fread("analyses/GRS_profiles/Wuttke2019_eGFR_PGSuncurated/profile.sscore.gz"),
  "Stroke_metaGRS"=fread("analyses/GRS_profiles/Stroke_metaGRS/profile.sscore.gz")
)
comp[PRS == "Khera2018_T2D_PGS000014", disease := "Type 2 Diabetes"]
comp[PRS == "Khera2018_CAD_PGS000013", disease := "Coronary Artery Disease"]
comp[PRS == "Khera2018_afib_PGS000016", disease := "Atrial Fibrillation"]
comp[PRS == "Mahajan2018_T2D_PGS000036", disease := "Type 2 Diabetes"]
comp[PRS == "Weng2017_afib_PGS000035", disease := "Atrial Fibrillation"]
comp[PRS == "Wuttke2019_eGFR_PGSuncurated", c("disease", "score_sum") := .("Chronic Kidney Disease", -1*score_sum)] # low eGFR = CKD
comp[PRS == "Stroke_metaGRS", disease := "Stroke"]
comp <- comp[, .(IID, PRS, disease, score_sum)]

# Standardise scores and adjust for 10 PCs
prs <- prs[pcs, on = .(IID)]
prs[, score := scale(score_sum), by=PRS]
prs[, score := lm(score ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +
                          PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals,
    by=PRS]
prs <- prs[,.(disease, PRS, IID, score)]

comp <- comp[pcs, on = .(IID)]
comp[, score := scale(score_sum), by=PRS]
comp[, score := lm(score ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +
                          PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals,
    by=PRS]
comp <- comp[,.(disease, PRS, IID, score)]

# Load associations for analysed PRSs
soma_assocs <- fread("analyses/pub/cardiometabolic_proteins/all_assocs.tsv")
soma_assocs <- soma_assocs[,.(PRS, Gene, Target, UniProt, Beta=Prot.Beta, 
                              L95=Prot.L95, U95=Prot.U95, Pvalue=Prot.Pvalue,
                              FDR=Prot.FDR)]
soma_assocs <- unique(soma_assocs)
soma_assocs[, FDR := FDR < 0.05]
soma_assocs[view, on = .(PRS=Display_name), c("PRS", "disease") := .(GRS_name, Display_name)]

# Load associations for PRSs we're comparing to:
rbindf <- function(...) { rbind(..., fill=TRUE) }
comp_assocs <- foreach(prs_id = unique(comp$PRS), .combine=rbindf) %do% {
  this_assocs <- fread(sprintf("analyses/univariate_associations/%s/somalogic_proteins/associations.tsv", prs_id))
  this_assocs[, c("PRS", "disease") := .(prs_id, comp[PRS == prs_id, disease][1])]
  this_assocs
}

# Aggregate associations to the protein level
soma_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
comp_assocs[soma_info, on = .(trait=variable), c("Gene", "Target", "UniProt") := 
            .(Gene.Name, TargetFullName, UniProt.Id.Current.at.Uniprot)]
comp_assocs <- comp_assocs[, .(Beta=mean(beta), L95=mean(l95), U95=mean(u95), Pvalue=mean(pval)),
                           by=.(PRS, disease, Gene, Target, UniProt)]

# Flip the associations for eGFR, because low eGFR == CKD.
comp_assocs[PRS == "Wuttke2019_eGFR_PGSuncurated", c("Beta", "L95", "U95") := .(-Beta, -U95, -L95)]

# Compare associations for each pair of comparable PRSs
comp_assocs <- merge(soma_assocs, comp_assocs, by=c("disease", "Gene", "Target", "UniProt"), suffixes=c("", ".comp"))

glist <- foreach(disease_id = unique(comp_assocs$disease), .combine=c) %do% {
  # Determine plot limits here:
  lim <- comp_assocs[disease == disease_id, c(pmin(min(L95.comp), min(L95)), pmax(max(U95.comp), max(U95)))]
  expand <- (lim[2] - lim[1]) * 0.05
  lim <- lim + c(-expand, expand)

	# Data.table for background ribbon
	rdt <- data.table(   x = c(lim[1], 0, lim[2]),
										ymin = c(lim[1], 0,      0),
										ymax = c(     0, 0, lim[2]))
  
  foreach(prs_id = comp_assocs[disease == disease_id, unique(PRS.comp)], .combine=c) %do% {
    # Make plot
    this_assocs <- comp_assocs[disease == disease_id & PRS.comp == prs_id]
    g <- ggplot(this_assocs) +
         aes(x=Beta, y=Beta.comp, xmin=L95, xmax=U95, ymin=L95.comp, ymax=U95.comp) +
         geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
				 geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
				 geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
				 geom_abline(intercept = 0, slope = 1, linetype=2, color="#737373") +
         geom_errorbarh(data=this_assocs[!(FDR)], alpha=0.2, height=0, size=0.3) +
         geom_errorbar(data=this_assocs[!(FDR)], alpha=0.2, width=0, size=0.3) +
         geom_point(data=this_assocs[!(FDR)], shape=19, alpha=0.5, size=0.3) +
         geom_errorbarh(data=this_assocs[(FDR)], colour="red", alpha=0.5, height=0, size=0.3) +
         geom_errorbar(data=this_assocs[(FDR)], colour="red", alpha=0.5, width=0, size=0.3) +
         geom_point(data=this_assocs[(FDR)], colour="red", shape=19, size=0.5) +
         scale_x_continuous(name=paste(disease_id, "PRS"), limits=lim, expand=c(0,0), breaks=pretty_breaks(n=3)) +
         scale_y_continuous(name=prs_id, limits=lim, expand=c(0,0), breaks=pretty_breaks(n=3)) +
         theme_bw() + 
				 theme(
					 axis.title=element_text(size=8), axis.text=element_text(size=8),
					 panel.grid=element_blank()
				 )

    return(list(g))
  }
}

assoc_plot <- plot_grid(glist[[1]], glist[[2]], glist[[3]], glist[[4]],
                        glist[[5]], glist[[6]], glist[[7]],
                        nrow=2, ncol=4)
ggsave(assoc_plot, width=7, height=3.2, file=sprintf("%s/pgs_sensitivity.png", out_dir), units="in", dpi=300) # raster image so it can be imported into inkscape

# Generate PDF version with the same dimensions without the points so we have pretty axis labels:
glist <- foreach(disease_id = unique(comp_assocs$disease), .combine=c) %do% {
  # Determine plot limits here:
  lim <- comp_assocs[disease == disease_id, c(pmin(min(L95.comp), min(L95)), pmax(max(U95.comp), max(U95)))]
  expand <- (lim[2] - lim[1]) * 0.05
  lim <- lim + c(-expand, expand)

  foreach(prs_id = comp_assocs[disease == disease_id, unique(PRS.comp)], .combine=c) %do% {
    # Make plot
    this_assocs <- comp_assocs[disease == disease_id & PRS.comp == prs_id]
    g <- ggplot(this_assocs) +
         aes(x=Beta, y=Beta.comp, xmin=L95, xmax=U95, ymin=L95.comp, ymax=U95.comp) +
         scale_x_continuous(name=paste(disease_id, "PRS"), limits=lim, expand=c(0,0), breaks=pretty_breaks(n=3)) +
         scale_y_continuous(name=prs_id, limits=lim, expand=c(0,0), breaks=pretty_breaks(n=3)) +
         theme_bw() + 
				 theme(
					 axis.title=element_text(size=8), axis.text=element_text(size=8),
					 panel.grid=element_blank()
				 )

    return(list(g))
  }
}

assoc_plot <- plot_grid(glist[[1]], glist[[2]], glist[[3]], glist[[4]],
                        glist[[5]], glist[[6]], glist[[7]],
                        nrow=2, ncol=4)
ggsave(assoc_plot, width=7, height=3.2, file=sprintf("%s/pgs_sensitivity.pdf", out_dir)) 

# Plot correlations between PGRS levels
comp <- merge(prs, comp, by=c("disease", "IID"), suffixes=c("", ".comp"))
glist <- foreach(disease_id = unique(comp_assocs$disease), .combine=c) %do% {
  # Determine plot limits here:
  xlim <- comp[disease == disease_id, max(abs(score))]
  expand <- xlim * 0.1
  xlim <- c(-xlim - expand, xlim + expand)

  ylim <- comp[disease == disease_id, c(min(score.comp), max(score.comp))]
  expand <- (ylim[2] - ylim[1]) * 0.05
  ylim <- ylim + c(-expand, expand)
  
  foreach(prs_id = comp[disease == disease_id, unique(PRS.comp)], .combine=c) %do% {
    # Make plot
    this_comp <- comp[disease == disease_id & PRS.comp == prs_id]
    g <- ggplot(this_comp, aes(x=score, y=score.comp)) +
				 geom_hline(yintercept=0, size=0.3, colour="#bdbdbd") +
				 geom_vline(xintercept=0, size=0.3, colour="#bdbdbd") +
				 geom_abline(intercept = 0, slope = 1, linetype=2, color="#bdbdbd", size=0.3) +
				 geom_point_rast(alpha=0.3, raster.width=3, raster.height=3, size=0.5) +
         annotate("text", color="red", x = -Inf, y = Inf, hjust = 0, vjust = 1,
				  				 label=paste0("r=", round(this_comp[,cor(score, score.comp)]*100)/100)) +
         scale_x_continuous(name=paste(disease_id, "PRS"), limits=xlim, expand=c(0,0), breaks=pretty_breaks(n=3)) +
         scale_y_continuous(name=prs_id, limits=ylim, expand=c(0,0), breaks=pretty_breaks(n=3)) +
				 theme_bw() +
				 theme(
           axis.text=element_text(size=8), axis.title=element_text(size=8),
				   panel.grid=element_blank()
         )

    return(list(g))
  }
}

comp_plot <- plot_grid(glist[[1]], glist[[2]], glist[[3]], glist[[4]],
                       glist[[5]], glist[[6]], glist[[7]],
                       nrow=2, ncol=4)
ggsave(comp_plot, width=7, height=3.2, file=sprintf("%s/pgs_correlations.pdf", out_dir))











