library(data.table)
library(pheatmap)
library(RColorBrewer)
library(foreach)
library(openxlsx)
library(ggplot2)
library(ggnewscale)
source("src/utilities/prot_pval.R")

out_dir <- "analyses/pub/cardiometabolic_proteins/review2"

# Load FDR significant associations
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])
prs_assocs <- prs_assocs[PRS.FDR < 0.05]

# Load mouse data
dat_dir <- sprintf("%s/mouse_data", out_dir)

ortho <- fread(sprintf("%s/mouse_orthologs.tsv", dat_dir))
hmdp_transcripts <- fread(sprintf("%s/hmdp_transcripts.tsv", dat_dir))
hmdp_cor <- fread(sprintf("%s/hmdp_trait_cor.tsv", dat_dir))
bxh_cor <- fread(sprintf("%s/bxh_trait_cor.tsv", dat_dir), colClasses=c("array_id"="character"))
qpcr <- as.data.table(read.xlsx("data/qPCR_mRNA_mice_AnnaCalkin.xlsx", sheet="Long format", startRow=3), check.names=TRUE)
time <- as.data.table(read.xlsx("data/qPCR_mRNA_mice_AnnaCalkin.xlsx", sheet="Longitudinal long format", startRow=3), check.names=TRUE)

# Build table about mouse orthologs
hmdp_transcripts[ortho, on = .(gene=hmdp.symbol), MGI.gene := mouse.symbol]
hmdp_transcripts <- hmdp_transcripts[,.(transcript_id=paste(unique(transcript), collapse=" ")), by=MGI.gene]
hmdp_transcripts <- hmdp_transcripts[transcript_id != ""]

bxh_transcripts <- unique(bxh_cor[,.(array_id, gene)])
bxh_transcripts[ortho, on = .(gene=bxh.symbol), MGI.gene := mouse.symbol]
bxh_transcripts <- bxh_transcripts[,.(transcript_id=paste(unique(array_id), collapse=" ")), by=MGI.gene]

ortho[hmdp_transcripts, on = .(mouse.symbol = MGI.gene), HMDP := transcript_id]
ortho[bxh_transcripts, on = .(mouse.symbol = MGI.gene), BXH := transcript_id]

ortho[mouse.symbol == "", mouse.symbol := "-"]
ortho[is.na(HMDP), HMDP := "-"]
ortho[is.na(BXH), BXH := "-"]

ortho <- ortho[order(mouse.symbol)][order(Gene)][order(PRS)]
fwrite(ortho[,.(PRS, UniProt, Gene, mouse.symbol, HMDP, BXH)], 
       sep="\t", quote=FALSE, file=sprintf("%s/mouse_orthologs.tsv", out_dir))

# For each gene to trait pair we want a single measure of the correlation. There are
# sometimes multiple transcripts per gene, and the datasets differ on the available 
# correlations. For the 'bxh_cor' correlations for all gene to trait pairs are 
# available (we computed them directly on another cluster where we have the data),
# but for 'hmdp_cor' correlations are only available where P < 0.05 and abs(Corr) > 0.1.
#
# In the case of the latter, to get a single correlation measure for each gene we must:
#
#   (1) Add back in "missing" correlation coefficients and P-values for gene to trait pairs
#       with abs(Corr) <= 0.1 and P >= 0.05. Since we don't know what these are, we impute 
#       them as Corr == 0 and P == 0.525, the midpoints for both the correlation and P-value 
#       within the cutoff ranges.
#
#   (2) Compute the average correlation coefficient and P-value across transcripts.
#
#   (3) Set gene to trait pairs with abs(Corr) <= 0.1 and P >= 0.05 to missing for consistency
#       with genes measured by a signle transcript.
#
# For the 'bxh_cor' we first set all correlation coefficients with abs(Corr) <= 0.01 or P > 0.05
# to missing so that we can repeat the process above and make the two heatmaps as comparable as 
# possible (i.e. to remove potential bias towards assocaitions in 'bxh_cor' and not 'hmdp_cor'
# introducted by us having all correlation coefficients).

# Make datasets consistent 
bxh_cor <- bxh_cor[abs(Corr) > 0.1 & P < 0.05]
bxh_cor[,  N := NULL]
setnames(hmdp_cor, c("tissue", "gene", "array_id", "trait", "Corr", "P")) 

# Add unifying information about MGI mouse gene symbol, protein, and list of
# transcripts for each gene.
hmdp_cor[ortho, on = .(gene=hmdp.symbol), 
         c("MGI.gene", "UniProt", "transcripts") := 
           .(i.mouse.symbol, i.UniProt, i.HMDP)]
bxh_cor[ortho, on = .(gene=bxh.symbol), 
         c("MGI.gene", "UniProt", "transcripts") := 
           .(i.mouse.symbol, i.UniProt, i.BXH)]

# How many transcripts per gene?
hmdp_cor[, N_transcripts := sapply(strsplit(transcripts, " "), length)]
bxh_cor[, N_transcripts := sapply(strsplit(transcripts, " "), length)]

# Fill in empty placeholder values
hmdp_empty <- hmdp_cor[,.N,by=.(tissue, UniProt, MGI.gene, trait, N_transcripts)][N != N_transcripts]
hmdp_empty <- hmdp_empty[, .(Corr = rep(0, N_transcripts - N), P = rep(0.5, N_transcripts - N)), by=.(tissue, UniProt, MGI.gene, trait)]
hmdp_cor <- rbind(hmdp_cor, hmdp_empty, fill=TRUE)

bxh_empty <- bxh_cor[,.N,by=.(tissue, UniProt, MGI.gene, trait, N_transcripts)][N != N_transcripts]
bxh_empty <- bxh_empty[, .(Corr = rep(0, N_transcripts - N), P = rep(0.5, N_transcripts - N)), by=.(tissue, UniProt, MGI.gene, trait)]
bxh_cor <- rbind(bxh_cor, bxh_empty, fill=TRUE)

# Average correlations across transcripts for each gene
hmdp_cor <- hmdp_cor[, .(Corr = mean(Corr), P = prot_pvalue(P, Corr)), by=.(tissue, UniProt, gene=MGI.gene, trait)]
bxh_cor <- bxh_cor[, .(Corr = mean(Corr), P = prot_pvalue(P, Corr)), by=.(tissue, UniProt, gene=MGI.gene, trait)]

# Drop genes where abs(Corr) <= 0.01 or P >= 0.05
hmdp_cor <- hmdp_cor[abs(Corr) > 0.01 & P < 0.05]
bxh_cor <- bxh_cor[abs(Corr) > 0.01 & P < 0.05]

# Filter each table to phenotypes of interest
hmdp_traits <- c( 
    # Lipids and cholesterol
    "HDL", "LDL_plus_VLDL", "TC", "TG",
    # Insulin resistance / sensitivity
    "Glucose", "Insulin", "glucose_to_insulin",
    # Body fat
    "BW", "NMR_BF_percentage"
  )

bxh_traits <- c( 
    # Atherosclerosis
    "Aortic Lesion", 
    # Lipids and cholesterol
    "HDL", "LDL+VLDL", "Total cholesterol", "Triglyceride", 
    # Insulin resistance / sensitivity
    "Glucose", "Insulin", "Glucose / Insulin",
    # Body fat
    "Leptin", "Weight", "100 x fat / weight"
  )

hmdp_cor <- hmdp_cor[trait %in% hmdp_traits]
bxh_cor <- bxh_cor[trait %in% bxh_traits]

# For the gene to trait correlations that are missing we want to add back in 0s
# so we can show these on the heatmaps, and distinguish from cases where the gene
# is not measured.
hmdp_map <- as.data.table(expand.grid(
  gene = ortho[HMDP != "-", unique(mouse.symbol)],
  trait = hmdp_traits, 
  tissue = unique(hmdp_cor$tissue)))
hmdp_map[ortho, on = .(gene=mouse.symbol), UniProt := i.UniProt]

hmdp_cor <- hmdp_cor[hmdp_map, on = .(tissue, trait, UniProt, gene)]
hmdp_cor[is.na(Corr), c("Corr", "P") := .(0, 0.5)]

bxh_map <- as.data.table(expand.grid(
  gene = ortho[BXH != "-", unique(mouse.symbol)],
  trait = bxh_traits, 
  tissue = unique(bxh_cor$tissue)))
bxh_map[ortho, on = .(gene=mouse.symbol), UniProt := i.UniProt]

bxh_cor <- bxh_cor[bxh_map, on = .(tissue, trait, UniProt, gene)]
bxh_cor[is.na(Corr), c("Corr", "P") := .(0, 0.5)]

# Add in NAs so the number of genes matches between datasets
hmdp_map <- as.data.table(expand.grid(
  gene = ortho[BXH != "-", unique(mouse.symbol)],
  trait = hmdp_traits, 
  tissue = unique(hmdp_cor$tissue)))
hmdp_map[ortho, on = .(gene=mouse.symbol), UniProt := i.UniProt]

hmdp_cor <- merge(hmdp_cor, hmdp_map, by = c("tissue", "trait", "UniProt", "gene"), all=TRUE)

bxh_map <- as.data.table(expand.grid(
  gene = ortho[HMDP != "-", unique(mouse.symbol)],
  trait = bxh_traits, 
  tissue = unique(bxh_cor$tissue)))
bxh_map[ortho, on = .(gene=mouse.symbol), UniProt := i.UniProt]

bxh_cor <- merge(bxh_cor, bxh_map, by = c("tissue", "trait", "UniProt", "gene"), all=TRUE)

# Combine tables and add PRS flag
all_cor <- rbind(HMDP=hmdp_cor, BXH=bxh_cor, idcol="dataset")
up_to_prs <- ortho[, .(PRS = paste(unique(PRS), collapse=";")), by=UniProt]
all_cor[up_to_prs, on = .(UniProt), PRS := i.PRS]

# Add in human gene
all_cor[ortho, on = .(gene=mouse.symbol), Human.Gene := i.Gene]

# Drop dataset label from tissues
all_cor[, tissue := gsub("HMDP_", "", tissue)]

# Output for paper
all_stat = copy(all_cor)
all_stat <- all_stat[, .(dataset, tissue, trait, Human.Gene, Mouse.gene=gene, Corr, P)]
fwrite(all_stat, file=sprintf("%s/mouse_trait_cor.csv", out_dir))

# Get left to right ordering of each protein, within each PRS group:
gene_order <- all_cor[,.(WD = sum(abs(Corr), na.rm=TRUE)/sum(!is.na(Corr))), by=.(PRS, UniProt, Human.Gene, gene)]
gene_order <- gene_order[order(-WD)][order(PRS)]

# Set up color palette
pal <- rev(brewer.pal(n = 8, name = "RdYlBu"))
pal <- c(pal[1:4], "#FFFFFF", pal[5:8])
pal <- colorRampPalette(pal)(255)

# Generate heatmaps:
for (dset in unique(all_cor$dataset)) {
  for (tiss in unique(all_cor[dataset == dset, tissue])) {
    bmat <- dcast(all_cor[dataset == dset & tissue == tiss], trait ~ gene, value.var="Corr")
    pmat <- dcast(all_cor[dataset == dset & tissue == tiss], trait ~ gene, value.var="P")

    bmat <- as.matrix(bmat, rownames=1)
    pmat <- as.matrix(pmat, rownames=1)

    # Order rows and columns
    if (dset == "HMDP") {
      bmat <- bmat[hmdp_traits, gene_order$gene]
      pmat <- pmat[hmdp_traits, gene_order$gene]
    } else {
      bmat <- bmat[bxh_traits, gene_order$gene]
      pmat <- pmat[bxh_traits, gene_order$gene]
    }

    # Arbitrary P-value signifier text
    ptext <- matrix("", nrow=nrow(pmat), ncol=ncol(pmat), dimnames=dimnames(pmat))
    ptext[pmat < 0.001] <- "*"

    # Create heatmap
    pheatmap(bmat, display_number=ptext, cluster_cols=FALSE, cluster_rows=FALSE,
             cellwidth=10, cellheight=7.2, fontsize_row=8, fontsize_col=8, fontsize_number=7,
             color=pal, breaks=seq(-1, 1, length=256), border_color=NA,
             file=sprintf("%s/%s_%s_heatmap.pdf", out_dir, dset, tiss))
  }
}

# Also reference heatmaps of PRS associations
prs_assocs <- prs_assocs[,. (Gene = strsplit(Gene, ",")[[1]]), by=.(PRS, UniProt, PRS.Beta)]
prs_assocs <- prs_assocs[gene_order$Human.Gene, on = .(Gene)]
prs_assocs <- prs_assocs[gene_order[,.(Human.Gene, Mouse.Gene=gene)], on = .(Gene=Human.Gene)]
prs_assocs <- unique(prs_assocs)
bmat <- dcast(prs_assocs, PRS ~ Mouse.Gene, value.var="PRS.Beta", fill=0)
bmat <- as.matrix(bmat, rownames=1)
bmat <- bmat[,gene_order$gene]
bmat[bmat < 0] <- -1
bmat[bmat > 0] <- 1

pheatmap(bmat, cluster_cols=FALSE, cluster_rows=FALSE,
         cellwidth=10, cellheight=7.2, fontsize_row=8, fontsize_col=8,
         color=c("#4393c3", "#ffffff", "#d6604d"), 
         border_color=NA, file=sprintf("%s/mouse_prs_ref_heatmap.pdf", out_dir))

# Generate boxplots for the dietary intervention qPCR data. We use a simple
# one-way anova to test, for each gene and tissue, whether the mean value
# is the same between the three diet groups (i.e. the null hypothesis is 
# mean(CHOW) == mean(HFD) == mean(WD)). These P-values are manually added
# to each plot in inkscape.
#
# According to the methods provided, the value in the qPCR column is the 
# delta delta Ct after standardisation to a housekeeping gene (Rplp0 
# expression in the tibialis anterior muscle (TA) and liver tissues, and
# cyclophilin A (Ppia) in the subcutaneous white adipose tissue (SC WAT)).
# This is apparently equivalent to the fold-change in gene expression 
# compared to the house keeping gene. 
#
# We log2-transform these fold-change values in the anova model in order 
# to account for the skewed non-normal nature of the data. On the fold-change
# scale a doubling of gene expression is a fold change of two, while 2-fold 
# decrease is a fold-change of 0.5. After log2 transform these values become
# 1 and -1 respectively

# First remove genes not analyses/no longer analysed here:
qpcr <- qpcr[Mouse.gene %in% ortho$mouse.symbol]
qpcr <- qpcr[Mouse.gene != "Shbg"]

# Remove genes with v. low expression in all samples/diets
qpcr <- qpcr[qPCR != 0] # 0 value reserved for "low expn" in all value cells.

# Remove failed samples
qpcr <- qpcr[!is.na(qPCR)]

# Write out for supp table
fwrite(qpcr, file=sprintf("%s/qpcr_raw.csv", out_dir))

# Get one-way anova statistics:
qpcr_aov <- foreach(gene = unique(qpcr$Mouse.gene), .combine=rbind) %do% {
  foreach(tissue = unique(qpcr$Tissue), .combine=rbind) %do% {
    dat <- qpcr[Tissue == tissue & Mouse.gene == gene]
    if (nrow(dat) > 0) {
      a1 <- summary(aov(log2(qPCR) ~ Diet, data=dat))
      data.table(Tissue = tissue, Human.gene = unique(dat$Human.gene), Mouse.gene = gene, 
                 Df = a1[[1]]["Diet", "Df"], SumSq = a1[[1]]["Diet", "Sum Sq"], MeanSq = a1[[1]]["Diet", "Mean Sq"],
                 Fstat = a1[[1]]["Diet", "F value"], P = a1[[1]][["Diet", "Pr(>F)"]])
    }
  }
}
fwrite(qpcr_aov, file=sprintf("%s/qpcr_aov.csv", out_dir))

# Get just p for ordering
qpcr_aov_p <- qpcr_aov[, .(Tissue, Human.gene, P)]

# Add PRS for ordering
g2p <- prs_assocs[, .(PRS = paste(unique(sort(PRS)), collapse=",")), by=Gene]
qpcr_aov_p[g2p, on = .(Human.gene=Gene), PRS := i.PRS]

# Add factor levels for within group ordering on plots
qpcr[, Diet := factor(Diet, levels=c("CHOW", "HFD", "WD"))]

# Duplicate the Diet column so we can use different fill values for the points
# and box plots (i.e. lighter colors in the boxplot).
qpcr[, Diet2 := Diet] 

# Generate box plots, ordering genes left to right by PRS, then by anova P-value.
g <- ggplot(qpcr[!is.na(qPCR) & Tissue == "Liver"],
            aes(x=factor(Human.gene, levels=qpcr_aov_p[Tissue == "Liver"][order(PRS, P), Human.gene]),
             y=qPCR)) +
  geom_hline(yintercept=1, color="#878787", linetype=2) +
  geom_boxplot(outlier.alpha=0, lwd=0.25, color="black", aes(fill=Diet)) +
  scale_fill_manual(name="Diet", values=c("CHOW"="#d9f0d3", "HFD"="#ffffbf", "WD"="#e7d4e8")) +
  new_scale_fill() +
  geom_point(position=position_jitterdodge(jitter.width=0), size=1, shape=21, color="black",
             aes(fill=Diet2)) +
  scale_fill_manual(name="Diet", values=c("CHOW"="#5aae61", "HFD"="#fee08b", "WD"="#762a83")) +
  #scale_y_log10(name="Fold-change", expand=expand_scale(mult=c(0.2, 0.2))) +
  scale_y_continuous(name="Fold-change", trans="log2", expand=expand_scale(mult=0.2),
                     breaks=c(0.25, 0.5, 1, 2, 4, 8)) +
  scale_x_discrete(name="", expand=expand_scale(add=0.5)) +
  theme_bw() +
  theme(axis.text=element_text(size=8, color="black"), axis.title=element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom")
ggsave(g, file=sprintf("%s/liver_qpcr.pdf", out_dir), units="in", width=4.4, height=2, useDingbats=FALSE)

g <- ggplot(qpcr[!is.na(qPCR) & Tissue == "SC WAT"],
            aes(x=factor(Human.gene, levels=qpcr_aov_p[Tissue == "SC WAT"][order(PRS, P), Human.gene]),
             y=qPCR)) +
  geom_hline(yintercept=1, color="#878787", linetype=2) +
  geom_boxplot(outlier.alpha=0, lwd=0.25, color="black", aes(fill=Diet)) +
  scale_fill_manual(name="Diet", values=c("CHOW"="#d9f0d3", "HFD"="#ffffbf", "WD"="#e7d4e8")) +
  new_scale_fill() +
  geom_point(position=position_jitterdodge(jitter.width=0), size=1, shape=21, color="black",
             aes(fill=Diet2)) +
  scale_fill_manual(name="Diet", values=c("CHOW"="#5aae61", "HFD"="#fee08b", "WD"="#762a83")) +
  scale_y_continuous(name="Fold-change", trans="log2", expand=expand_scale(mult=0.2),
                                          breaks=c(0.25, 0.5, 1, 2, 4)) +
  scale_x_discrete(name="", expand=expand_scale(add=0.5)) +
  theme_bw() +
  theme(axis.text=element_text(size=8, color="black"), axis.title=element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom")
ggsave(g, file=sprintf("%s/sc_wat_qpcr.pdf", out_dir), units="in", width=4.4, height=2, useDingbats=FALSE)

g <- ggplot(qpcr[!is.na(qPCR) & Tissue == "TA" & Time.point == "12 weeks"],
            aes(x=factor(Human.gene, levels=qpcr_aov_p[Tissue == "TA"][order(PRS, P), Human.gene]),
             y=qPCR)) +
  geom_hline(yintercept=1, color="#878787", linetype=2) +
  geom_boxplot(outlier.alpha=0, lwd=0.25, color="black", aes(fill=Diet)) +
  scale_fill_manual(name="Diet", values=c("CHOW"="#d9f0d3", "HFD"="#ffffbf", "WD"="#e7d4e8")) +
  new_scale_fill() +
  geom_point(position=position_jitterdodge(jitter.width=0), size=1, shape=21, color="black",
             aes(fill=Diet2)) +
  scale_fill_manual(name="Diet", values=c("CHOW"="#5aae61", "HFD"="#fee08b", "WD"="#762a83")) +
  scale_y_continuous(name="Fold-change", trans="log2", expand=expand_scale(mult=0.2),
                                          breaks=c(0.25, 0.5, 1, 2, 4, 8)) +
  scale_x_discrete(name="", expand=expand_scale(add=0.5)) +
  theme_bw() +
  theme(axis.text=element_text(size=8, color="black"), axis.title=element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom")
ggsave(g, file=sprintf("%s/ta_qpcr.pdf", out_dir), units="in", width=2.8, height=2, useDingbats=FALSE)

# CRYZL1 expression was additionally measured at dietary intervetions for 4 weeks, 12 weeks, and 16 weeks
# (all data shown above is for a 12 week dietary intervention). We will plot boxplots over time, showing
# linear model fits across time for each diet, as well as boxplots at each time point.

# Is the mean expression of each diet at each time point different for CRYZL1?
aov_p <- function(...) {  summary(aov(...))[[1]][["Pr(>F)"]][1] }
cryzl1_aov_p <- time[, .(P=aov_p(log2(qPCR) ~ Diet)), by=Time.point]

# Is there are a dose response effect over time?
time[, week := as.numeric(gsub(" .*", "", Time.point))]
cryzl1_lm <- time[, .(
                    Intercept=coef(summary(lm(log2(qPCR) ~ week)))[1,1],
                    Beta=coef(summary(lm(log2(qPCR) ~ week)))[2,1],
                    P=coef(summary(lm(log2(qPCR) ~ week)))[2,4]
                  ), by=Diet]

# Build data table for plotting lm fit. We need to do 2^(Beta) and 2^(Intercept)
# to get the values on the fold-change scale, the same as the rest of the data,
# then need to manually calculate the value at weeks = 4, 12, and 16 because
# week is a factor, not a continuous value, for the box plots to work.
# We plot these as points, so that a line can be later added through them in 
# inkscape, because ggplot really doesnt want to make this plot work (I can't seem
# to get boxplots by a discrete value on with continuous x).
lmfit <- copy(cryzl1_lm)
lmfit <- rbind("4"=lmfit, "12"=lmfit, "16"=lmfit, idcol="week")
lmfit[, log2y := Intercept + Beta * as.numeric(week)]
lmfit[, y := 2^log2y] 
lmfit[, week := factor(week, levels=c(4,8,12,16))]

# Duplicate the Diet column so we can use different fill values for the points
# and box plots (i.e. lighter colors in the boxplot).
time[, Diet2 := Diet] 

g <- ggplot(time[!is.na(qPCR)], 
            aes(x=factor(week, levels=c(4,8,12,16)), y=qPCR)) +
  geom_hline(yintercept=1, color="#878787", linetype=2) +
  geom_boxplot(outlier.alpha=0, lwd=0.25, color="black", aes(fill=Diet)) +
  scale_fill_manual(name="Diet", values=c("CHOW"="#d9f0d3", "HFD"="#ffffbf", "WD"="#e7d4e8")) +
  new_scale_fill() +
  geom_point(position=position_jitterdodge(jitter.width=0), size=1, shape=21, color="black",
             aes(fill=Diet2)) +
  scale_fill_manual(name="Diet", values=c("CHOW"="#5aae61", "HFD"="#fee08b", "WD"="#762a83")) +
  geom_point(data=lmfit, inherit.aes=FALSE, aes(x=week, y=y, color=Diet), 
             position=position_jitterdodge(jitter.width=0),
             shape=19, size=0.1) +
  scale_color_discrete(guide=FALSE) +
  scale_y_continuous(name="Fold-change", trans="log2", expand=expand_scale(mult=0.2),
                     breaks=c(0.5,1,2,4,8,16,32)) +
  scale_x_discrete(name="", drop=FALSE) +
  theme_bw() +
  theme(axis.text=element_text(size=8, color="black"), axis.title=element_text(size=10),
        axis.title.x=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom")
ggsave(g, file=sprintf("%s/cryzl1_times_qpcr.pdf", out_dir), units="in", width=2.8, height=2, useDingbats=FALSE)

