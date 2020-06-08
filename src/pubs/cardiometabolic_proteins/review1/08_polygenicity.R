library(data.table)
library(foreach)
library(ggplot2)
library(ggrastr)
library(scales)
library(cowplot)

source("src/utilities/prot_pval.R")

out_dir="analyses/pub/cardiometabolic_proteins/review1"

soma_assocs <- fread("analyses/pub/cardiometabolic_proteins/all_assocs.tsv")
soma_assocs <- soma_assocs[Prot.FDR < 0.05, .(PRS, Gene, Target, UniProt, Beta=Prot.Beta, L95=Prot.L95, U95=Prot.U95, Pvalue=Prot.Pvalue)]
soma_assocs[PRS == "Coronary Artery Disease", PRS := "CAD_metaGRS"]
soma_assocs[PRS == "Chronic Kidney Disease", PRS := "CKD_2019"]
soma_assocs[PRS == "Type 2 Diabetes", PRS := "T2D_2018"]
soma_assocs <- unique(soma_assocs)

soma_info = fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
pqtl_removed <- foreach(idx = soma_assocs[,.I], .combine=rbind) %do% {
  prs_id = soma_assocs[idx, PRS]
  gene_id = soma_assocs[idx, Gene]
  apt_assocs <- foreach(var_id = soma_info[Gene.Name == gene_id, variable], .combine=rbind) %do% {
    fread(sprintf("analyses/grs_pqtl_removed/%s/%s/associations.tsv", prs_id, var_id))
  }
  apt_assocs[, .(PRS=prs_id, Gene=gene_id, Beta=mean(beta), L95=mean(l95), U95=mean(u95), Pvalue=prot_pvalue(pval, beta))]
}

has_pqtls <- fread("analyses/sensitivity_associations/qtl_prob_dosage_adjusted/CAD_metaGRS/somalogic_proteins/qtl_grs_associations.tsv")
has_pqtls <- has_pqtls[variable != "grs"]
has_pqtls <- unique(has_pqtls[,.(trait)])
has_pqtls[soma_info, on = .(trait = variable), Gene := Gene.Name]
has_pqtls <- unique(has_pqtls[,.(Gene)])
soma_assocs[, has_pQTL := FALSE]
soma_assocs[has_pqtls, on = .(Gene), has_pQTL := TRUE]

# First, generate new plot which condenses what was previously in Figure S4 into a single panel:
ggdat <- merge(soma_assocs, pqtl_removed, by=c("PRS", "Gene"), suffixes=c("", ".nopqtl"))
ggdat[PRS == "CAD_metaGRS", PRS := "Coronary Artery Disease"]
ggdat[PRS == "CKD_2019", PRS := "Chronic Kidney Disease"]
ggdat[PRS == "T2D_2018", PRS := "Type 2 Diabetes"]

# data.table for ribbon
rdt <- ggdat[, .(x=c(min(c(L95, L95.nopqtl)), 0,
                    max(c(U95, U95.nopqtl)))*1.05,
                ymax=c(min(c(L95, L95.nopqtl)), 0,
                       max(c(U95, U95.nopqtl)))*1.05,
                ymin=c(0,0,0)),
             by=PRS]

g <- ggplot(ggdat, aes(x = Beta, y = Beta.nopqtl, xmin = L95, ymin = L95.nopqtl, xmax = U95, ymax = U95.nopqtl)) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
  geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#737373") +
  geom_errorbarh(data=ggdat[!(has_pQTL)], height=0, alpha=0.5, size=0.5, color="#636363") +
  geom_errorbar(data=ggdat[!(has_pQTL)], width=0, alpha=0.5, size=0.5, color="#636363") +
  geom_point(data=ggdat[!(has_pQTL)], shape = 19, size=1, color="#636363") +
  geom_errorbarh(data=ggdat[(has_pQTL)], height=0, alpha=0.5, size=0.5, color="#1f78b4") +
  geom_errorbar(data=ggdat[(has_pQTL)], width=0, alpha=0.5, size=0.5, color="#1f78b4") +
  geom_point(data=ggdat[(has_pQTL)], shape = 19, size=1, color="#1f78b4") +
  facet_wrap(~ PRS, ncol=3, scales="free") +
  scale_x_continuous(name = "Beta (95% CI) for the PRS", expand=c(0,0)) +
  scale_y_continuous(name = "Beta (95% CI)\nexcluding pQTL regions", expand=c(0,0)) +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    strip.text=element_text(size=10), panel.grid=element_blank(),
  )
ggsave(g, width=7.2, height=2.4, file=sprintf("%s/pqtl_removal_sensitivity.pdf", out_dir), useDingbats=FALSE)

# Load in PRS variants
prs <- fread("analyses/pub/cardiometabolic_proteins/review1/filtered_score_files.txt")

# Load in chunked association results
cis_windows <- fread(sprintf("%s/cis_regions.txt", out_dir))

chunk_assocs <- fread(sprintf("%s/one_chunk_assocs.txt", out_dir))
polygenicity <- fread(sprintf("%s/progressive_chunk_leaveout_assocs.txt", out_dir))
gwas <- fread(sprintf("%s/sig_prot_score_variant_pQTL_effects.txt", out_dir))

# Define chunks
gwas[, chunk_10MB_num := floor(pos/10e6)]
prs[, chunk_10MB_num := floor(pos/10e6)]

# Plot distribution of N score variants per chunk
prs_n_per_chunk <- prs[,.N,by=.(PRS, chr, chunk_10MB_num)]
g <- ggplot(prs_n_per_chunk, aes(x=PRS, y=N)) +
  geom_boxplot() +
  scale_y_log10(name="Score variants per 10MB chunk") +
  xlab("") +
  theme_bw() +
  theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8)
  )
ggsave(g, file=sprintf("%s/prs_variants_per_chunk.pdf", out_dir))

# Define cumulative numeric positions for manhattan plots
cumul_pos <- gwas[,.(min=min(pos), max=max(pos)),by=chr]
cumul_pos[,pos_offset := c(0, cumsum(as.numeric(max)+1)[-22])]
cumul_pos[, chr_label := pos_offset + max/2]

gwas[cumul_pos, on = .(chr), cumul_pos := pos + pos_offset]

cis_windows[cumul_pos, on = .(chr), cumulative_start := start + pos_offset]
cis_windows[cumul_pos, on = .(chr), cumulative_end := end + pos_offset]
cis_windows[cumul_pos, on = .(chr), cumulative_TSS := TSS + pos_offset]

chunk_assocs[, chunk_mid := chunk_start + (chunk_end - chunk_start)/2]
chunk_assocs[cumul_pos, on = .(chr), chunk_mid_cumul := chunk_mid + pos_offset]

# table to define cis-window ribbons for plotting
ggcis <- melt(cis_windows, id.vars=c("Gene"), measure.vars=c("cumulative_start", "cumulative_end"), value.name="cumul_pos")
ggcis[, ymin := -Inf]
ggcis[, ymax := Inf]

# Ribbons for full score effects
assoc_ribbon <- soma_assocs[, .(PRS, Gene, L95, U95)]
assoc_ribbon <- assoc_ribbon[, .(x=c(0, cumul_pos[22, max + pos_offset]), ymin=L95, ymax=U95), by=.(PRS, Gene)]
assoc_ribbon[polygenicity[,.(max=max(n_chunks)),by=PRS], on = .(PRS), n_chunks := max]
assoc_ribbon[x == 0, n_chunks := 1]

pqtl_removed_ribbon <- pqtl_removed[, .(PRS, Gene, L95, U95)]
pqtl_removed_ribbon <- pqtl_removed_ribbon[, .(x=c(0, cumul_pos[22, max + pos_offset]), ymin=L95, ymax=U95), by=.(PRS, Gene)]
pqtl_removed_ribbon[polygenicity[,.(max=max(n_chunks)),by=PRS], on = .(PRS), n_chunks := max]
pqtl_removed_ribbon[x == 0, n_chunks := 1]

# Function to make full polygenicity plot for a PRS to protein pair
plot_polygenicity <- function(prs_id, gene_id) {
  prs_name <- gsub("_.*", "", prs_id)

  # Manhattan plot of pQTLs associations with the protein for variants in the PRS
	pqtl_manhattan <- ggplot(gwas[Gene == gene_id][prs[PRS == prs_id], on = .(chr, pos)]) +
		aes(x=cumul_pos, y=-log10(pval), color=factor(chr %% 2)) +
    geom_vline(xintercept=cis_windows[Gene == gene_id, cumulative_TSS], color="#fd8d3c") +
		geom_point(shape=19) +
		geom_hline(yintercept=-log10(5e-8), linetype=2, color="red") +
		scale_colour_manual(guide=FALSE, values=c("0"="#9ecae1", "1"="#3182bd")) +
		scale_x_continuous(name="", breaks=cumul_pos$chr_label, labels=gsub("21", "", cumul_pos$chr), expand=c(0,0)) +
		scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0, 0.05))) +
		ggtitle(sprintf("pQTL effects for %s PRS variants on plasma %s levels", prs_name, gene_id)) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(),
			axis.title = element_text(size=10), axis.text=element_text(size=8), title=element_text(size=8)
		)

	# Plot similar manhattan plots for each 10 MB score region
	chunk_manhattan <- ggplot(chunk_assocs[PRS == prs_id & Gene == gene_id]) +
    aes(x=chunk_mid_cumul, y=-log10(Pvalue), color=factor(chr %% 2)) +
    geom_vline(xintercept=cis_windows[Gene == gene_id, cumulative_TSS], color="#fd8d3c") +
		geom_point(shape=19) +
		geom_hline(yintercept=soma_assocs[PRS == prs_id & Gene == gene_id, -log10(Pvalue)], linetype=2, color="#7a0177") + 
		geom_hline(yintercept=pqtl_removed[PRS == prs_id & Gene == gene_id, -log10(Pvalue)], linetype=2, color="#c51b8a") + 
		scale_colour_manual(guide=FALSE, values=c("0"="#9ecae1", "1"="#3182bd")) +
		scale_x_continuous(name="", breaks=cumul_pos$chr_label, labels=gsub("21", "", cumul_pos$chr), expand=c(0,0)) +
		scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0, 0.05))) +
		ggtitle(sprintf("Association of each 10MB chunk of %s PRS with plasma %s levels", prs_name, gene_id)) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(),
			axis.title = element_text(size=10), axis.text=element_text(size=8), title=element_text(size=8)
		)

	# Plot beta effects for each 10 MB score region
	beta_chunks <- ggplot(chunk_assocs[PRS == prs_id & Gene == gene_id]) +
		aes(x=chunk_mid_cumul, y=Beta, ymin=L95, ymax=U95, color=factor(chr %% 2)) +
    geom_vline(xintercept=cis_windows[Gene == gene_id, cumulative_TSS], color="#fd8d3c") +
		geom_ribbon(data=assoc_ribbon[PRS == prs_id & Gene == gene_id], inherit.aes=FALSE, aes(x=x, ymin=ymin, ymax=ymax), fill="#7a0177", alpha=0.2) +
		geom_ribbon(data=pqtl_removed_ribbon[PRS == prs_id & Gene == gene_id], inherit.aes=FALSE, aes(x=x, ymin=ymin, ymax=ymax), fill="#c51b8a", alpha=0.2) +
		geom_hline(yintercept=soma_assocs[PRS == prs_id & Gene == gene_id, Beta], linetype=2, color="#7a0177") + 
		geom_hline(yintercept=pqtl_removed[PRS == prs_id & Gene == gene_id, Beta], linetype=2, color="#c51b8a") + 
		geom_errorbar(width=0, size=0.3, alpha=0.5) +
		geom_point(shape=19) +
		geom_hline(yintercept=0, linetype=2, color="#636363") +
		scale_colour_manual(guide=FALSE, values=c("0"="#3182bd", "1"="#08519c")) +
		scale_x_continuous(name="", breaks=cumul_pos$chr_label, labels=gsub("21", "", cumul_pos$chr), expand=c(0,0)) +
		scale_y_continuous(name="Beta (95% CI)") +
		ggtitle(sprintf("SD change in plasma %s levels per SD increase in each 10MB chunk of %s PRS", gene_id, prs_name)) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(),
			axis.title = element_text(size=10), axis.text=element_text(size=8), title=element_text(size=8)
		)

	# Plot change in association P-value as number of chunks is decreased
	pval_change <- ggplot(polygenicity[PRS == prs_id & Gene == gene_id]) + 
    aes(x=n_chunks, y=-log10(Pvalue)) +
		geom_line(colour="#3182bd") +
		geom_hline(yintercept=soma_assocs[PRS == prs_id & Gene == gene_id, -log10(Pvalue)], linetype=2, color="#7a0177") +
		geom_hline(yintercept=pqtl_removed[PRS == prs_id & Gene == gene_id, -log10(Pvalue)], linetype=2, color="#c51b8a") +
		scale_x_reverse(name="Number of 10MB chunks included in score", expand=expansion(mult=c(0.01,0.01))) +
		scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0.01, 0.05))) +
		ggtitle(sprintf("%s PRS to plasma %s association with successive removal of 10MB chunks by P-value", prs_name, gene_id)) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(),
			axis.title = element_text(size=10), axis.text=element_text(size=8), title=element_text(size=8),
			axis.title.x = element_text(size=8, colour="#636363")
		)
  
  # Plot change in association beta as number of chunks is decreased
  beta_change <- ggplot(polygenicity[PRS == prs_id & Gene == gene_id]) +
		aes(x=n_chunks, y=Beta, ymin=L95, ymax=U95) + 
		geom_ribbon(data=assoc_ribbon[PRS == prs_id & Gene == gene_id], inherit.aes=FALSE, aes(x=n_chunks, ymin=ymin, ymax=ymax), fill="#7a0177", alpha=0.2) +
		geom_ribbon(data=pqtl_removed_ribbon[PRS == prs_id & Gene == gene_id], inherit.aes=FALSE, aes(x=n_chunks, ymin=ymin, ymax=ymax), fill="#c51b8a", alpha=0.2) +
		geom_hline(yintercept=soma_assocs[PRS == prs_id & Gene == gene_id, Beta], linetype=2, color="#7a0177") + 
		geom_hline(yintercept=pqtl_removed[PRS == prs_id & Gene == gene_id, Beta], linetype=2, color="#c51b8a") + 
		geom_ribbon(fill="#3182bd", alpha=0.5) +
		geom_line(colour="#08519c") +
		geom_hline(yintercept=0, linetype=2, color="#636363") +
		scale_x_reverse(name="Number of 10MB chunks included in score", expand=expansion(add=c(1,1))) +
		scale_y_continuous(name="Beta (95% CI)") +
		ggtitle(sprintf("SD change in plasma %s levels per SD increase of %s PRS with successive chunk removal", gene_id, prs_name)) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(),
			axis.title = element_text(size=10), axis.text=element_text(size=8), title=element_text(size=8),
			axis.title.x = element_text(size=8, colour="#636363")
		)

	# Aggregate plots
	g <- plot_grid(ncol=1,
		pqtl_manhattan, 
		chunk_manhattan,
		beta_chunks,
		pval_change,
		beta_change
	)
	ggsave(g, width=7.2, height=7.2, units="in", dpi=100, file=sprintf("%s/%s_to_%s_polygenicity.png", out_dir, prs_id, gene_id))
}

foreach(idx = soma_assocs[,.I]) %do% {
  gc()
  plot_polygenicity(soma_assocs[idx, PRS], soma_assocs[idx, Gene])
}

# Create nice versions for the T2D PRS to SHBG association for the paper
# Load in cis-pQTLs - so that we can obtain the cis-pvalue threshold for SHBG
cis_hier <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
cis_hier[soma_info, on = .(SOMAMER_ID), c("Gene", "Target", "UniProt") := .(Gene.Name, TargetFullName, UniProt.Id.Current.at.Uniprot)]

# Manhattan plot of pQTLs associations with the protein for variants in the PRS
pqtl_manhattan <- ggplot(gwas[Gene == "SHBG"][prs[PRS == "T2D_2018"], on = .(chr, pos), nomatch=0]) +
  aes(x=cumul_pos, y=-log10(pval), color=factor(chr %% 2)) +
  geom_vline(xintercept=cis_windows[Gene == "SHBG", cumulative_TSS], color="#fd8d3c") +
  geom_point_rast(shape=19, raster.width=10, raster.height=1, size=1, raster.dpi=300) +
  geom_hline(yintercept=-log10(1.5e-11), linetype=2, color="red") + # Trans pQTL significance threshold from Sun et al. 2018
  #geom_hline(yintercept=-log10(cis_hier[Gene == "SHBG", max(P)]), linetype=2, color="#fd8d3c") + # Cis pQTL significance threshold.
  scale_colour_manual(guide=FALSE, values=c("0"="#9ecae1", "1"="#3182bd")) +
  scale_x_continuous(name="Chromosome", breaks=cumul_pos$chr_label, labels=gsub("21", "", cumul_pos$chr), expand=c(0,0.0001)) +
  scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0, 0.1))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), 
    panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
    axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
  )
ggsave(pqtl_manhattan, width=7.2, height=1, file=sprintf("%s/T2D_PRS_to_SHBG_pqtls.pdf", out_dir))

# Plot similar manhattan plots for each 10 MB score region
chunk_manhattan <- ggplot(chunk_assocs[PRS == "T2D_2018" & Gene == "SHBG"]) +
  aes(x=chunk_mid_cumul, y=-log10(Pvalue), color=factor(chr %% 2)) +
  geom_vline(xintercept=cis_windows[Gene == "SHBG", cumulative_TSS], color="#fd8d3c") +
  geom_point_rast(shape=19, raster.width=12, raster.height=1, size=2, raster.dpi=300) +
  geom_hline(yintercept=soma_assocs[PRS == "T2D_2018" & Gene == "SHBG", -log10(Pvalue)], linetype=2, color="#7a0177") + 
  geom_hline(yintercept=-log10(0.05), linetype=2, color="red") + 
  scale_colour_manual(guide=FALSE, values=c("0"="#9ecae1", "1"="#3182bd")) +
  scale_x_continuous(name="Chromosome", breaks=cumul_pos$chr_label, labels=gsub("21", "", cumul_pos$chr), expand=c(0,0)) +
  scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0, 0.1))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), 
    panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
    axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
  )
ggsave(chunk_manhattan, width=7.2, height=1, file=sprintf("%s/T2D_PRS_to_SHBG_10MB.pdf", out_dir))

# Plot change in association P-value as number of chunks is decreased
chunks_per_prs <- polygenicity[,.(total_chunks=max(n_chunks)), by=PRS]
polygenicity[chunks_per_prs, on = .(PRS), prop_removed := (total_chunks - n_chunks)/total_chunks]
pval_change <- ggplot(polygenicity[PRS == "T2D_2018" & Gene == "SHBG"]) + 
  aes(x=prop_removed*100, y=-log10(Pvalue)) +
  geom_line(colour="#3182bd") +
  geom_hline(yintercept=-log10(0.05), linetype=2, color="red") +
  scale_x_continuous(name="Proportion of genome removed from T2D PRS", limits=c(0, 100), expand=expansion(mult=c(0.01,0.01))) +
  scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0.02, 0.05))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), 
    panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
    axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
  )
ggsave(pval_change, width=7.2, height=1, file=sprintf("%s/T2D_PRS_to_SHBG_polygenicity.pdf", out_dir))

# How many chunks do we need to remove to attenuate association?
polygenicity[chunks_per_prs, on = .(PRS), n_removed := total_chunks - n_chunks]
attenuation <- rbind(idcol="boundary",
  "upper"=polygenicity[Pvalue < 0.05, .(N=max(n_removed)+1), by=.(PRS, Gene)],
  "lower"=polygenicity[Pvalue > 0.05, .(N=min(n_removed)), by=.(PRS, Gene)])
attenuation[chunks_per_prs, on = .(PRS), pct := N/total_chunks]

# Plotting labels and order
attenuation[, factor_order := seq_len(.N), by=boundary]
attenuation[PRS == "T2D_2018", display_name := "Type 2 Diabetes"]
attenuation[PRS == "CAD_metaGRS", display_name := "Coronary Artery\nDisease"]
attenuation[PRS == "CKD_2019", display_name := "Chronic Kidney\nDisease"]
attenuation[,Gene.Label := Gene]
attenuation[Gene == "C5orf38", Gene.Label := "CEI"]
attenuation[Gene == "PDE4D", Gene.Label := "PDE4A/D*"]

# Get median for plotting
median <- attenuation[,.(median=median(N)),by=.(boundary, display_name)]

g <- ggplot(attenuation, aes(x=factor(factor_order), y=N, fill=boundary)) +
  geom_hline(data=median, aes(yintercept=median, colour=boundary), linetype=2) +
  geom_col(position="identity") +
  scale_colour_manual(name="", guide=FALSE, values=c("lower"="#7a0177", "upper"="#ae017e")) +
  scale_fill_manual(name="", values=c("lower"="#7a0177", "upper"="#ae017e"), labels=c(
    "lower"="Smallest N where PRS to protein association P > 0.05",
    "upper"="Largest N where PRS to protein association P < 0.05"
  )) +
  scale_x_discrete(name="", labels=attenuation[boundary=="lower", structure(Gene.Label, names=factor_order)]) +
  scale_y_continuous(name="Number of 10MB chunks", breaks=c(0, 30, 60, 90, 120), 
                     labels=c("0", "30 (10%)", "60 (21%)", "90 (31%)", "120 (41%)")) +
  facet_grid(. ~ display_name, scales="free_x", space="free_x") +
  theme_bw() +
  theme(
		axis.title = element_text(size=8), axis.text=element_text(size=8),
    axis.text.x=element_text(size=8, angle = 90, hjust = 1, vjust=0.5),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.position="bottom", legend.direction="horizontal",
    legend.text=element_text(size=8)
	)
ggsave(g, width=7.2, height=2.8, file=sprintf("%s/chunk_attenuation_complex.pdf", out_dir))


# Simpler version for paper figure
g <- ggplot(attenuation[boundary == "lower"], aes(x=factor(factor_order), y=pct*100)) +
  geom_col(position="identity", fill="#7a0177") +
  scale_x_discrete(name="", labels=attenuation[boundary=="lower", structure(Gene.Label, names=factor_order)]) +
  scale_y_continuous(name="Percentage of genome") +
  facet_grid(. ~ display_name, scales="free_x", space="free_x") +
  theme_bw() +
  theme(
    axis.title = element_text(size=8), axis.text=element_text(size=8),
    axis.text.x=element_text(size=8, angle = 90, hjust = 1, vjust=0.5),
    strip.text=element_text(size=8),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
  )
ggsave(g, width=7.2, height=1.8, file=sprintf("%s/chunk_attenuation_simple.pdf", out_dir))

# Show also the association between CKD PRS and FTMT

# Manhattan plot of pQTLs associations with the protein for variants in the PRS
pqtl_manhattan <- ggplot(gwas[Gene == "FTMT"][prs[PRS == "CKD_2019"], on = .(chr, pos), nomatch=0]) +
  aes(x=cumul_pos, y=-log10(pval), color=factor(chr %% 2)) +
  geom_vline(xintercept=cis_windows[Gene == "FTMT", cumulative_TSS], color="#fd8d3c") +
  geom_point_rast(shape=19, raster.width=10, raster.height=1, size=1, raster.dpi=300) +
  geom_hline(yintercept=-log10(1.5e-11), linetype=2, color="red") + # Trans pQTL significance threshold from Sun et al. 2019
  #geom_hline(yintercept=-log10(cis_hier[Gene == "FTMT", max(P)]), linetype=2, color="#fd8d3c") + # Cis pQTL significance threshold.
  scale_colour_manual(guide=FALSE, values=c("0"="#9ecae1", "1"="#3182bd")) +
  scale_x_continuous(name="Chromosome", breaks=cumul_pos$chr_label, labels=gsub("21", "", cumul_pos$chr), expand=c(0,0.0001)) +
  scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0, 0.1))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), 
    panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
    axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
  )
ggsave(pqtl_manhattan, width=7.2, height=1, file=sprintf("%s/CKD_PRS_to_FTMT_pqtls.pdf", out_dir))

# Plot similar manhattan plots for each 10 MB score region
chunk_manhattan <- ggplot(chunk_assocs[PRS == "CKD_2019" & Gene == "FTMT"]) +
  aes(x=chunk_mid_cumul, y=-log10(Pvalue), color=factor(chr %% 2)) +
  geom_vline(xintercept=cis_windows[Gene == "FTMT", cumulative_TSS], color="#fd8d3c") +
  geom_point_rast(shape=19, raster.width=12, raster.height=1, size=2, raster.dpi=300) +
  geom_hline(yintercept=soma_assocs[PRS == "CKD_2019" & Gene == "FTMT", -log10(Pvalue)], linetype=2, color="#7a0177") + 
  geom_hline(yintercept=-log10(0.05), linetype=2, color="red") + 
  scale_colour_manual(guide=FALSE, values=c("0"="#9ecae1", "1"="#3182bd")) +
  scale_x_continuous(name="Chromosome", breaks=cumul_pos$chr_label, labels=gsub("21", "", cumul_pos$chr), expand=c(0,0)) +
  scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0, 0.1))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), 
    panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
    axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
  )
ggsave(chunk_manhattan, width=7.2, height=1, file=sprintf("%s/CKD_PRS_to_FTMT_10MB.pdf", out_dir))

# Plot change in association P-value as number of chunks is decreased
chunks_per_prs <- polygenicity[,.(total_chunks=max(n_chunks)), by=PRS]
polygenicity[chunks_per_prs, on = .(PRS), prop_removed := (total_chunks - n_chunks)/total_chunks]
pval_change <- ggplot(polygenicity[PRS == "CKD_2019" & Gene == "FTMT"]) + 
  aes(x=prop_removed*100, y=-log10(Pvalue)) +
  geom_line(colour="#3182bd") +
  geom_hline(yintercept=-log10(0.05), linetype=2, color="red") +
  scale_x_continuous(name="Proportion of genome removed from CKD PRS", limits=c(0, 100), expand=expansion(mult=c(0.01,0.01))) +
  scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0.02, 0.05))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), 
    panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
    axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
  )
ggsave(pval_change, width=7.2, height=1, file=sprintf("%s/CKD_PRS_to_FTMT_polygenicity.pdf", out_dir))
