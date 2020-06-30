library(data.table)
library(foreach)
library(ggplot2)
library(ggrepel)
library(ggrastr)

out_dir <- "analyses/pub/cardiometabolic_proteins/review1"

# Load view file
grs_info <- fread("views/cardiometabolic.txt")

# Load in associations and filter to those with FDR < 0.05
soma_assocs <- fread("analyses/pub/cardiometabolic_proteins/all_assocs.tsv")
soma_assocs <- soma_assocs[Prot.FDR < 0.05]

# Load NGS assocs
ngs_assocs <- foreach(grs_idx = grs_info$GRS_name, .combine=rbind) %do% {
  dt <- fread(sprintf("analyses/univariate_associations/%s/olink_ngs/associations.tsv", grs_idx))
  dt[, PRS := grs_idx]
  dt
}
ngs_assocs[grs_info, on = .(PRS=GRS_name), PRS := Display_name]

# get mapping between Somalogic and Olink proteins
soma_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
soma_info <- soma_info[,.(Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, Gene=Gene.Name)]
soma_info <- unique(soma_info)
soma_info <- soma_info[unique(soma_assocs[,.(Target, UniProt, Gene)]), on = .(Target, UniProt, Gene)]

ngs_info <- fread("analyses/processed_traits/olink_ngs/trait_info.tsv")
ngs_info <- ngs_info[,.(UniProt, Gene, variable)]
ngs_info <- ngs_info[,.(UniProt = strsplit(UniProt, "_")[[1]]), by=.(Gene, variable)]
soma_info[ngs_info, on = .(UniProt), Olink := variable]

# Filter both association results to proteins measured on both platfomrs
shared <- soma_info[!is.na(Olink)]

ngs_assocs <- ngs_assocs[shared, on = .(trait=Olink), nomatch=0]
setnames(ngs_assocs, "trait", "Olink")

soma_assocs <- soma_assocs[shared, on = .(Target, UniProt, Gene), nomatch=0]

ngs_assocs <- ngs_assocs[soma_assocs[,.(PRS, Gene)], on = .(PRS, Gene)]

# build table comparing associations
soma_assocs <- unique(soma_assocs[,.(PRS, Target, UniProt, Gene, beta=Prot.Beta, l95=Prot.L95, u95=Prot.U95, pval=Prot.Pvalue, fdr=Prot.FDR)])
ngs_assocs <- ngs_assocs[,.(PRS, Target, UniProt, Gene, beta, l95, u95, pval)]

comp <- merge(soma_assocs, ngs_assocs, by=c("PRS", "Target", "UniProt", "Gene"), suffixes=c(".soma", ".olink"))

# Compare PRS to protein associations
col_map <- c(
  "Coronary Artery Disease"="#7570b3",
  "Type 2 Diabetes"="#d95f02",
  "Chronic Kidney Disease"="#1b9e77"
)

plot_min <- comp[, min(c(l95.soma, l95.olink))]
plot_max <- comp[, max(c(u95.soma, u95.olink))]
expand <- (plot_max - plot_min) * 0.1
plot_min <- plot_min - expand
plot_max <- plot_max + expand

rdt <- data.table(   x = c(plot_min, 0, 0, 0, plot_max),
                  ymin = c(0, 0, 0, plot_min, plot_min),
                  ymax = c(plot_max, plot_max, 0, 0, 0))

cor_assocs <- comp[, cor.test(beta.soma, beta.olink)]
cor_label <- sprintf("Pearson r = %s\n[95%% CI: %s-%s]",
                     round(cor_assocs$estimate*100)/100,
                     round(cor_assocs$conf.int[1]*100)/100,
                     round(cor_assocs$conf.int[2]*100)/100)

g1 <- ggplot(comp) +
  aes(x = beta.soma, y = beta.olink,
      xmin = l95.soma, ymin = l95.olink,
      xmax = u95.soma, ymax = u95.olink,
      label = Gene,
      color=PRS, fill=PRS) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
  geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#737373") +
  geom_errorbarh(height=0, alpha=0.5, size=0.5) +
  geom_errorbar(width=0, alpha=0.5, size=0.5) +
  geom_point(shape = 21, size=1.3, color="#00000000") +
  geom_text_repel(color="black", size=2, nudge_x=0.01) +
  geom_text(inherit.aes=FALSE, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1, label = cor_label) +
  scale_x_continuous(name = "Beta (95% CI) somalogic platform", expand=c(0,0)) +
  scale_y_continuous(name = "Beta (95% CI) olink NGS platform", expand=c(0,0)) +
  scale_color_manual(name = "PRS", values=col_map,
                     guide = guide_legend(direction="horizontal", nrow=3)) +
  scale_fill_manual(name = "PRS", values=col_map,
                    guide = guide_legend(direction="horizontal", nrow=3)) +
  coord_fixed(ratio=1) +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g1, width=3.6, height=3.6, file=sprintf("%s/soma_olink_NGS_association_compare.pdf", out_dir), useDingbats=FALSE)

# As forest plots
comp <- rbind(idcol="platform", "SomaLogic"=soma_assocs, "Olink NGS"=ngs_assocs, fill=TRUE)
comp[,Gene := factor(Gene, levels = comp[platform == "SomaLogic"][order(-beta), Gene])] # order genes by effect size
comp[PRS == "Chronic Kidney Disease", PRS := "Chronic Kidney\nDisease"]
comp[PRS == "Coronary Artery Disease", PRS := "Coronary Artery\nDisease"]
comp[,platform := factor(platform, levels=c("SomaLogic", "Olink NGS"))]

g2 <- ggplot(comp, aes(x=Gene, y=beta, ymin=l95, ymax=u95, colour=platform, fill=platform)) +
  geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_errorbar(width=0, alpha=0.5, size=0.5, position=position_dodge(width=0.8)) +
  geom_point(shape = 21, size=1.3, color="#00000000", position=position_dodge(width=0.8)) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Beta (95% CI)") +
  scale_colour_manual(name="Platform", values=c("SomaLogic"="#e08214", "Olink NGS"="#8073ac")) +
  scale_fill_manual(name="Platform", values=c("SomaLogic"="#e08214", "Olink NGS"="#8073ac")) +
  facet_grid(. ~ PRS, space="free_x", drop=TRUE, shrink=TRUE, scales="free_x") +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid=element_blank(), legend.position="bottom", legend.title=element_text(size=10),
    legend.text=element_text(size=8), strip.text=element_text(size=8, face="bold"), 
    strip.background=element_blank()
  )
ggsave(g2, width=7.2, height=3, file=sprintf("%s/soma_olink_NGS_association_forest.pdf", out_dir), useDingbats=FALSE)


# Load olink protein levels measured by PCR and NGS to compare directly
ngs <- fread("analyses/processed_traits/olink_ngs/traits.tsv")
olink <- fread("analyses/processed_traits/olink_proteins/traits.tsv")
olink_info <- fread("analyses/processed_traits/olink_proteins/trait_info.tsv")

# Identify shared proteins, then filter values
shared2 <- olink_info[shared, on = .(UniProt), nomatch=0][panel != "neu"]

ngs <- ngs[shared2[,.(Target, Gene, UniProt, Olink)], on = .(variable=Olink)]
olink <- olink[shared2[,.(Target, Gene, UniProt, variable)], on = .(variable)]

comp <- merge(ngs, olink, by=c("IID", "Target", "Gene", "UniProt"), suffixes=c(".ngs", ".pcr"))


# plot comparison of protein levels
cor_label <- comp[, .(
  label = sprintf("Pearson r = %s\n[95%% CI: %s-%s]",
                  round(cor.test(value.pcr, value.ngs)$estimate*100)/100,
                  round(cor.test(value.pcr, value.ngs)$conf.int[1]*100)/100,
                  round(cor.test(value.pcr, value.ngs)$conf.int[2]*100)/100)
), by=Gene]

g3 <- ggplot(comp, aes(x=value.pcr, y=value.ngs)) +
  geom_hline(yintercept=0, size=0.3, colour="#bdbdbd") +
  geom_vline(xintercept=0, size=0.3, colour="#bdbdbd") +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#bdbdbd", size=0.3) +
  geom_point_rast(alpha=0.3, raster.width=3, raster.height=3, size=1) +
  geom_text(data=cor_label, color="red", size=2, x = -Inf, y = Inf, hjust = 0, vjust = 1, aes(label=label)) +
  facet_wrap(~ Gene, ncol=2) +
  scale_x_continuous(name="Olink PCR levels (old)", breaks=c(-3, 0, 3), limits=c(-4, 4)) +
  scale_y_continuous(name="Olink NGS levels (new)", breaks=c(-3, 0, 3), limit=c(-4,4)) +
  theme_bw() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),
        strip.text=element_text(size=10, face="bold"),
        strip.background=element_blank(), panel.grid=element_blank())
ggsave(g3, width=3.6, height=4, file=sprintf("%s/olink_seqtype_compare.pdf", out_dir))








