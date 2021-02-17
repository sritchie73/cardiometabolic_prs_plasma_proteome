library(data.table)
library(ggplot2)
library(scales)

# Load all polygenicity results
polygenicity <- fread("analyses/pub/cardiometabolic_proteins/review2/polygenicity.txt")
one_chunk_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/one_chunk_assocs.txt")

# Extract the LD blocks contributing to polygenicity of each association
one_chunk_assocs <- one_chunk_assocs[order(P)][order(Gene)][order(PRS)]
one_chunk_assocs[, block_num := 1:.N, by=.(PRS, Gene)]
one_chunk_assocs <- one_chunk_assocs[polygenicity[,.(PRS, Gene, LD_blocks_removed)], on = .(PRS, Gene, block_num <= LD_blocks_removed)]

# Build table of chromosome position sizes - we will plot these first so chromosomes
# are the right size (and so we don't have to plot a tonne of empty LD blocks)
chrom_blocks <- one_chunk_assocs[,.(block_start=min(block_start), block_end=max(block_end)), by=chr]
chrom_blocks <- polygenicity[,.(chr=1:22), by=.(PRS, Gene)][chrom_blocks, on = .(chr)]

# Create y-axis labels and order
chrom_blocks[, ytext := sprintf("%s PGS to %s", gsub("_PRS", "", PRS), Gene)]
chrom_blocks[, ytext := gsub("C5orf38", "CEI", ytext)]
chrom_blocks[, ytext := gsub("PDE4D", "PDE4D/A", ytext)]
chrom_blocks <- chrom_blocks[polygenicity[order(-pct_removed), .(PRS, Gene)], on = .(PRS, Gene)]
chrom_blocks <- rbind(chrom_blocks[PRS == "T2D_PRS"], chrom_blocks[PRS == "CAD_PRS"], chrom_blocks[PRS == "IS_PRS"], chrom_blocks[PRS == "CKD_PRS"])
chrom_blocks[, ytext := factor(ytext, levels=unique(ytext))]

one_chunk_assocs[chrom_blocks, on = .(PRS, Gene), ytext := i.ytext]

# Get cis location
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
cis_loc <- unique(prs_assocs[,.(Gene, chr, start)])
cis_loc <- cis_loc[unique(chrom_blocks[,.(PRS, Gene, ytext)]), on = .(Gene)]
cis_loc <- rbind(cis_loc, data.table(Gene = "PDE4A", chr=19, start=10527449, PRS="CKD_PRS", ytext="CKD PGS to PDE4D/A"))
cis_loc[, start := as.integer(start)]
cis_loc <- cis_loc[chr != "X"]
cis_loc[,chr := as.integer(chr)]

# Make plot
g <- ggplot(chrom_blocks) +
  aes(xmin = block_start, xmax = block_end, fill=-log10(P)) +
  geom_rect(data=one_chunk_assocs[-log10(P) <= 5], ymin = 0, ymax = 1, color=NA, size=0) +
  geom_rect(data=one_chunk_assocs[-log10(P) > 5], ymin = 0, ymax = 1, color=NA, fill="#4a1486", size=0) +
  geom_point(data=cis_loc, inherit.aes=FALSE, aes(x=start), y=1, color="black", shape=18, size=1.2) + 
  geom_point(data=cis_loc, inherit.aes=FALSE, aes(x=start), y=0, color="black", shape=18, size=1.2) + 
  geom_rect(ymin = 0, ymax = 1, fill=NA, color="black", size=0.2) +
  scale_x_continuous(name="Chromosome", expand=c(0,0)) +
  scale_y_continuous(name="", expand=c(0,0)) +
  scale_fill_gradient(low="#fed976", high="#e31a1c") +
  facet_grid(ytext ~ chr, scales="free_x", space="free", switch="both") +
  theme_bw() + theme(
    strip.text.x.bottom=element_text(size=5), strip.text.y.left=element_text(size=6, angle=0, vjust=0.5, hjust=1),
    strip.switch.pad.grid=unit(0.2, "mm"), strip.background=element_blank(),
    panel.spacing.x = unit(0.5, "mm"), panel.spacing.y = unit(0.5, "mm"), 
    panel.border=element_blank(), panel.background=element_blank(), panel.grid=element_blank(),
    axis.text.x=element_blank(), axis.title.x=element_text(size=7), axis.ticks=element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"),
    legend.position="bottom", legend.title=element_text(size=7), legend.text=element_text(size=6)
  )
ggsave(g, width=18.3, height=18.3, units="cm", file="analyses/pub/cardiometabolic_proteins/review3/polygenicity.pdf")



