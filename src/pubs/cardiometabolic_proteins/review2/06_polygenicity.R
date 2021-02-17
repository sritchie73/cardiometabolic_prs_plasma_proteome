library(data.table)
library(foreach)
library(ggplot2)
library(ggrastr)
library(scales)
library(cowplot)
library(doMC)

source("src/utilities/prot_pval.R")
source("src/utilities/flip_strand.R")

# Warning: don't request a full node here because if you do for some reason
# none of the parallelism works - recommend requesting 20 cores on skylake-himem
# and parallelising foreach loops across 10 cores.
registerDoMC(10)
setDTthreads(20)

# Load HES data
source("src/pubs/cardiometabolic_proteins/review2/load_HES.R")

# Filter to primary+secondary endpoints
hes <- hes[code_type == "primary+secondary"]

# Load list of cardiometabolic events for prevalent event exclusion
source("src/pubs/cardiometabolic_proteins/review2/cardiometabolic_events.R")

# Exclude people with prevalent events
prev <- unique(hes[prevalent == 1 & phenotype %in% cardiometabolic, .(IID)])
hes <- hes[!prev, on = .(IID)]

# Load SomaLogic aptamer levels
soma <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")

# Adjust for batch
batch <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")
soma <- soma[batch, on = .(IID), nomatch=0]
soma[, value := lm(value ~ factor(batch))$residuals, by=variable]
soma[, soma_ivt_adj_batch := scale(value)]
soma <- soma[,.(IID, variable, soma_ivt_adj_batch)]

# Drop prevalent events
soma <- soma[!prev, on = .(IID)]

# Drop withdrawn sample
idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
idmap <- fread(idmap_file)
soma <- soma[IID %in% idmap$Affymetrix_gwasQC_bl]

# Load protein information.
sinfo <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

# Filter to aptamers passing QC
sinfo <- sinfo[Type == "Protein"]

# Select columns
sinfo <- sinfo[Type == "Protein", .(variable, SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot,
                                    Gene=Gene.Name, Entrez_id=Entrez.Gene.ID, chr, start, end)]

# Fix bad entries (Aptamers for the same target with different/missing gene/uniprot information)
sinfo[Target == "14-3-3 protein family", UniProt := "P61981|Q04917"]
sinfo[Target == "Induced myeloid leukemia cell differentiation protein Mcl-1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("Q07820", "MCL1", "4170", "1", "150547027", "150552214")]
sinfo[Target == "Protein delta homolog 1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("P80370", "DLK1", "8788", "14", "101193202", "101201467")]
sinfo[Target == "Stromal cell-derived factor 1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("P48061", "CXCL12", "6387", "10", "44865601", "44880545")]

# Get prs-associated proteins and aptamers
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- prs_assocs[FDR < 0.05]

# Load pcs
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

#################################################
# Get variants in PRS
#################################################

if (file.exists("analyses/pub/cardiometabolic_proteins/review2/filtered_score_files.txt")) {
  scores <- fread("analyses/pub/cardiometabolic_proteins/review2/filtered_score_files.txt")
} else {

	# Load PRS and filter to set of variants used (ultimately just need the chr and pos)
	scores <- rbind(idcol="PRS",
		CAD_PRS=fread("data/GRS_resources/CAD_metaGRS/grs_weights.txt.gz"),
		IS_PRS=fread("data/GRS_resources/Stroke_metaGRS/grs_weights.txt.gz"),
		T2D_PRS=fread("data/GRS_resources/T2D_2018/grs_weights.txt.gz"),
		CKD_PRS=fread("data/GRS_resources/CKD_2019/grs_weights.txt.gz")
	)

	scores <- scores[effect_allele != flip_strand(other_allele)]

	# Load in bim files to find matches
	scores[, match := FALSE]
	for (chrIdx in 1:22) {
		gc()
		bim <- fread(sprintf("analyses/processed_genotypes/impute_chr%s_interval_filtered.bim", chrIdx), header=FALSE)
		setnames(bim, c("chr", "rsid", "cm", "pos", "A1", "A2"))

		scores[bim, on = .(chr, pos, effect_allele=A1, other_allele=A2), match := TRUE]
		scores[bim, on = .(chr, pos, effect_allele=A2, other_allele=A1), match := TRUE]

		scores[bim, on = .(chr, pos), c("effect_allele", "other_allele") :=
					 .(ifelse(match, effect_allele, flip_strand(effect_allele)),
						 ifelse(match, other_allele, flip_strand(other_allele)))]
		scores[bim, on = .(chr, pos, effect_allele=A1, other_allele=A2), match := TRUE]
		scores[bim, on = .(chr, pos, effect_allele=A2, other_allele=A1), match := TRUE]
		rm(bim)
	}

	# Drop variants without any match:
	scores <- scores[(match)]

	# Drop any multi-allelic score variants
	mult <- scores[,.N,by=.(PRS, pos)][N > 1]
	scores <- scores[!mult, on = .(PRS, pos)]

	# Write out
	scores[, match := NULL]
	fwrite(scores, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/filtered_score_files.txt")
}

#########################################################
# Get pQTL summary statistics for all tests of interest
#########################################################


if (file.exists("analyses/pub/cardiometabolic_proteins/review2/pQTLs_for_prs_variants.txt")) {
  prs_pqtls <- fread("analyses/pub/cardiometabolic_proteins/review2/pQTLs_for_prs_variants.txt")
  tests <- prs_assocs[!unique(prs_pqtls[,.(PRS, Target)]), on = .(PRS, Target)]
} else {
  tests <- prs_assocs
}

if (nrow(tests) > 0) {
	prs_apt_qtls <- foreach(test_idx = tests[,.I], .combine=rbind) %:% 
		foreach(this_chr = 1:22, .combine=rbind) %dopar% {
			this_prs <- prs_assocs[test_idx, PRS]
			this_target <- prs_assocs[test_idx, Target]
			this_uniprot <- prs_assocs[test_idx, UniProt]
			this_gene <- prs_assocs[test_idx, Gene]
			this_seqid <- prs_assocs[test_idx, Aptamer]
			this_aptid <- sinfo[Aptamer == this_seqid, SOMAMER_ID]

			this_prs_vars <- scores[PRS == this_prs & chr == this_chr, .(pos)]
			this_apt_qtls <- fread(sprintf("data/full_pQTL_summary_stats/%s/%s_chrom_%s_meta_1.tbl.gz", this_aptid, this_aptid, this_chr))
			this_apt_qtls <- this_apt_qtls[this_prs_vars, on = .(position=pos), nomatch=0]

			this_apt_qtls <- this_apt_qtls[, .(position, Effect, P=10^`log(P)`)]
			
			iter_info <- data.table(Target=this_target, UniProt=this_uniprot, Gene=this_gene, Aptamer=this_seqid, PRS=this_prs, chr=this_chr)
			cbind(iter_info, this_apt_qtls)
	}

	prs_apt_qtls <- prs_apt_qtls[, .(P=prot_pvalue(P, Effect)), by = .(Target, UniProt, Gene, PRS, chr, position)]

  # Add flag for cis regions
  cis_windows <- sinfo[, .(chr = strsplit(chr, "\\|")[[1]], TSS=strsplit(start, "\\|")[[1]]), by=.(Target, UniProt, Gene, Aptamer)]
  cis_windows <- unique(cis_windows[, .(Target, UniProt, Gene, chr, TSS)])
  cis_windows <- cis_windows[, .(chr = strsplit(chr, ";")[[1]], TSS=strsplit(TSS, ";")[[1]]), by=.(Target, UniProt, Gene)]
  cis_windows <- cis_windows[chr != "X" & chr != "Y"]
  cis_windows[, chr := as.integer(chr)]
  cis_windows[, TSS := as.integer(TSS)]
	cis_windows[, start := pmax(TSS - 1e6, 0)]
	cis_windows[, end := TSS + 1e6]

	complex_ld <- data.table(
		region_chr=c(5, 6, 8, 11),
		region_start=c(44000000, 25000000, 8000000, 45000000),
		region_end=c(51500000, 33500000, 12000000, 57000000),
		region_name=c("r1", "MHC", "r3", "r4")
	)

	cis_windows[complex_ld, on = .(chr = region_chr, start > region_start, end < region_end),
						c("start", "end") := .(region_start, region_end)]
	cis_windows[complex_ld, on = .(chr = region_chr, start < region_end, end > region_end),
						c("start") := .(region_start)]
	cis_windows[complex_ld, on = .(chr = region_chr, end > region_start, start < region_start),
						c("end") := .(region_end)]

  prs_apt_qtls[, cis := FALSE]
  prs_apt_qtls[cis_windows, on = .(Target, UniProt, Gene, chr, position >= start, position <= end), cis := TRUE]


  if (exists("prs_pqtls")) {
    prs_pqtls <- rbind(prs_pqtls, prs_apt_qtls)
  } else {
    prs_pqtls <- prs_apt_qtls
  }

	fwrite(prs_pqtls, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/pQTLs_for_prs_variants.txt")
}

##########################################################
# Test each PRS chunk for association with protein levels
##########################################################

# Load chunked PRS
chunked_scores <- fread("analyses/pub/cardiometabolic_proteins/review2/chunked_scores.sscore.gz")

if (!("PRS" %in% names(chunked_scores))) {
  chunked_scores[, PRS := gsub("_PRS.*", "_PRS", score_chunk)]  
  chunked_scores[, ld_block := as.integer(gsub(".*_", "", score_chunk))]
  chunked_scores[, score_chunk := NULL]
  chunked_scores <- chunked_scores[!is.na(ld_block)] # 3 variants in AF PRS on different chromosomes outside defined LD blocks.
  chunked_scores <- chunked_scores[,.(IID, PRS, ld_block, score_sum)]

  fwrite(chunked_scores, sep="\t", quote=FALSE, compress="gzip", file="analyses/pub/cardiometabolic_proteins/review2/chunked_scores.sscore.gz")
}

# Filter to individauls without prevalent events and with somalogic data
chunked_scores <- chunked_scores[IID %in% unique(soma$IID)]

# Test each chunk for association with the aptamer of interest
tests <- prs_assocs[,.(PRS, Target, UniProt, Gene, Aptamer)]
tests[sinfo, on = .(Aptamer), aptvar := variable]

apt_one_chunk_assocs <- foreach(this_test = tests[,.I], .combine=rbind) %do% {
  this_prs <- tests[this_test, PRS]
  this_target <- tests[this_test, Target]
  this_gene <- tests[this_test, Gene]
  this_uniprot <- tests[this_test, UniProt]
  this_aptamer <- tests[this_test, Aptamer]
  this_aptvar <- tests[this_test, aptvar]

  foreach(this_chunk = chunked_scores[PRS == this_prs, unique(ld_block)], .combine=rbind) %dopar% {
    dat <- chunked_scores[PRS == this_prs & ld_block == this_chunk]
    dat <- dat[pcs, on = .(IID), nomatch=0]
    dat[, pc_adj_score := scale(lm(score_sum ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +
                                                 PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals)]
    dat <- dat[soma[variable == this_aptvar], on = .(IID), nomatch=0]

    l1 <- lm(soma_ivt_adj_batch ~ pc_adj_score, data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)

    data.table(PRS = this_prs, ld_block = this_chunk, Target = this_target, UniProt = this_uniprot, 
               Gene = this_gene, Aptamer = this_aptamer, Beta = cf[2,1], SE = cf[2,2], L95 = ci[2, 1],
               U95 = ci[2,2], P = cf[2,4])
  }
}

# Average across multiple aptamers for the given protein
one_chunk_assocs <- apt_one_chunk_assocs[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)), 
                                         by=.(PRS, ld_block, Target, UniProt, Gene)]

# Load information about each LD block
ld_blocks <- fread("data/Berisa_etal_2016/EUR_1000G_ind_ld_blocks.bed")
ld_blocks[, chr := as.integer(gsub("chr", "", chr))]
ld_blocks[, block_size := stop - start]
ld_blocks[, ld_block := .I]
ld_blocks <- ld_blocks[,.(ld_block, chr, block_start = start, block_end = stop, block_size)]

one_chunk_assocs <- merge(ld_blocks, one_chunk_assocs, by="ld_block")

fwrite(one_chunk_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/one_chunk_assocs.txt")

# Also write out combined version with aptamer information
comb <- merge(one_chunk_assocs, apt_one_chunk_assocs, by = c("PRS", "ld_block", "Target", "UniProt", "Gene"), suffixes=c("", ".aptamer"))
fwrite(comb, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/one_chunk_assocs_with_aptamer.txt")

##########################################################
# Progressive leave one out analysis
##########################################################

one_chunk_assocs <- one_chunk_assocs[order(-abs(Beta))][order(P)][order(Target)][order(PRS)]

tests <- prs_assocs[,.(PRS, Target, UniProt, Gene, Aptamer)]
tests[sinfo, on = .(Aptamer), aptvar := variable]
apt_leaveout_assocs <- foreach(this_test = tests[,.I], .combine=rbind) %dopar% {
  this_prs <- tests[this_test, PRS]
  this_target <- tests[this_test, Target]
  this_gene <- tests[this_test, Gene]
  this_uniprot <- tests[this_test, UniProt]
  this_aptamer <- tests[this_test, Aptamer]
  this_aptvar <- tests[this_test, aptvar]

  foreach(chunk_idx = one_chunk_assocs[PRS == this_prs & Target == this_target, .I], .combine=rbind) %do% {
    leaveout <- one_chunk_assocs[PRS == this_prs & Target == this_target][1:chunk_idx, ld_block]
    bp_removed <- one_chunk_assocs[PRS == this_prs & Target == this_target][1:chunk_idx, sum(block_size)]
    pct_removed <- bp_removed / ld_blocks[,sum(block_size)]

    dat <- chunked_scores[PRS == this_prs & !(ld_block %in% leaveout)]
    if (nrow(dat) == 0) {
      return(NULL)
    }
    dat <- dat[, .(score_sum=sum(score_sum)), by=.(IID)]

    dat <- dat[pcs, on = .(IID), nomatch=0]
    dat[, pc_adj_score := scale(lm(score_sum ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +
                                                 PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals)]
    dat <- dat[soma[variable == this_aptvar], on = .(IID), nomatch=0]

    l1 <- lm(soma_ivt_adj_batch ~ pc_adj_score, data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)

    data.table(PRS = this_prs, LD_blocks_removed=chunk_idx, BP_removed=bp_removed, pct_removed=pct_removed,
               Target = this_target, UniProt = this_uniprot, 
               Gene = this_gene, Aptamer = this_aptamer, Beta = cf[2,1], SE = cf[2,2], L95 = ci[2, 1],
               U95 = ci[2,2], P = cf[2,4])
  }
}

# Average across multiple aptamers for the given protein
leaveout_assocs <- apt_leaveout_assocs[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)), 
                                          by=.(PRS, LD_blocks_removed, BP_removed, pct_removed, Target, UniProt, Gene)]

leaveout_assocs <- rbind(
  unique(prs_assocs[,.(PRS, LD_blocks_removed=0, BP_removed=0, pct_removed=0, Target, UniProt, Gene, Beta, SE, L95, U95, P)]),
  leaveout_assocs
)

fwrite(leaveout_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/leaveout_chunk_assocs.txt")

# Also write out combined version with aptamer information
comb <- merge(leaveout_assocs, apt_leaveout_assocs, by = c("PRS", "Target", "UniProt", "Gene", "LD_blocks_removed", "BP_removed", "pct_removed"), suffixes=c("", ".aptamer"), all.x=TRUE)
fwrite(comb, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/leaveout_assocs_with_aptamer.txt")

##########################################################
# Compute polygenicity
##########################################################

leaveout_assocs <- leaveout_assocs[order(LD_blocks_removed)][order(Target)][order(PRS)]

polygenicity <- leaveout_assocs[P > 0.05, .SD[1], by=.(PRS, Target, UniProt, Gene)]
polygenicity <- polygenicity[, .(PRS, Target, UniProt, Gene, LD_blocks_removed, BP_removed, pct_removed)]

fwrite(polygenicity, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/polygenicity.txt")

##########################################################
# Plot polygenicity
##########################################################

# Order proteins in decreasing order of association strength
prot_order <- unique(prs_assocs[, .(PRS, Gene, Beta, P)])
prot_order[, PRS := factor(PRS, levels=c("T2D_PRS", "CAD_PRS", "IS_PRS", "CKD_PRS"))]
prot_order <- prot_order[order(-abs(Beta))][order(P)][order(PRS)]
prot_order[, prot_order := .I]
polygenicity[prot_order, on = .(PRS, Gene), prot_order := i.prot_order]
polygenicity <- polygenicity[order(prot_order)]

g <- ggplot(polygenicity) +
  aes(x=factor(prot_order), y=pct_removed*100) +
  geom_col(position="identity", fill="#7a0177") +
  scale_x_discrete(name="", breaks=polygenicity$prot_order, labels=polygenicity$Gene) +
  scale_y_continuous(name="Percentage of genome") +
  facet_grid(. ~ PRS, scales="free_x", space="free_x") +
  theme_bw() +
  theme(
    axis.title = element_text(size=8), axis.text=element_text(size=8),
    axis.text.x=element_text(size=8, angle = 90, hjust = 1, vjust=0.5),
    strip.text=element_text(size=8),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
  )
ggsave(g, width=7.2, height=1.8, file="analyses/pub/cardiometabolic_proteins/review2/polygenicity_barchart.pdf")

##########################################################
# Plot pQTLs vs chunk assocs vs polygenicity
##########################################################

# Get cis pQTL threshold for each gene
cis_thresh <- fread("analyses/mendelian_randomisation/pqtls/aptamer_cis_hierarchical_pthresh.tsv")
cis_thresh[sinfo, on = .(SOMAMER_ID), c("Target", "UniProt", "Gene") := .(Target, UniProt, Gene)]
cis_thresh <- cis_thresh[, .(threshold=prot_pvalue(threshold, rep(1, .N))), by = .(Target, UniProt, Gene)]

# Define cumulative genome position
last_pos <- ld_blocks[,.SD[which.max(ld_block)], by=chr]
cumul_pos <- foreach(idx = 1:last_pos[,.N], .combine=c) %do% {
  1 + last_pos[1:idx, sum(block_end)]
}
cumul_pos <- c(0, cumul_pos)

one_chunk_assocs <- one_chunk_assocs[order(ld_block)][order(Gene)][order(PRS)]
tests <- unique(prs_assocs[,.(PRS, Target, UniProt, Gene, chr, TSS=start, P)])
foreach(this_test = tests[,.I]) %dopar% {
  this_prs <- tests[this_test, PRS]
  this_gene <- tests[this_test, Gene]
  this_uniprot <- tests[this_test, UniProt]
  this_chr <- tests[this_test, chr]
  this_tss <- as.integer(tests[this_test, TSS])

  # Plot overlap of pQTLs and PRS variants
  this_pqtls <- prs_pqtls[PRS == this_prs & Gene == this_gene]
  this_pqtls <- this_pqtls[order(position)][order(chr)]
  this_pqtls[, cumul_pos := .I]
  cumul_pos_tss <- this_pqtls[chr == this_chr][which.min(abs(position - this_tss)), cumul_pos]
  xticks <- this_pqtls[, .(chr_lab_pos=median(cumul_pos)), by=chr]
  this_cis_thresh <- cis_thresh[Gene == this_gene, threshold]
	pqtl_manhattan <- ggplot(this_pqtls) +
		aes(x=cumul_pos, y=-log10(P), color=factor(chr %% 2)) +
		geom_vline(xintercept=cumul_pos_tss, color="#fd8d3c") +
		geom_point_rast(shape=19, raster.width=10, raster.height=1, size=1, raster.dpi=300) +
		geom_hline(yintercept=-log10(1.5e-11), linetype=2, color="red") + # Trans pQTL significance threshold from Sun et al. 2018
		geom_hline(yintercept=-log10(this_cis_thresh), linetype=2, color="#fd8d3c") + # Cis pQTL significance threshold.
		scale_colour_manual(guide=FALSE, values=c("0"="#9ecae1", "1"="#3182bd")) +
		scale_x_continuous(name="Chromosome", breaks=xticks$chr_lab_pos, labels=gsub("21", "", 1:22), expand=c(0,0.0001)) +
		scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0.05, 0.1))) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(),
			panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
			axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
		)

	# Plot associations for each independent LD block
  this_one_chunk_assocs <- one_chunk_assocs[PRS == this_prs & Gene == this_gene]
  this_one_chunk_assocs[, xpos := cumul_pos[chr] + block_start + (block_size)/2, by=ld_block] 
  this_prs_assoc_p <- tests[this_test, P]
  xticks <- (cumul_pos[1:22] + cumul_pos[1:22 + 1])/2
  cumul_pos_tss <- cumul_pos[as.integer(this_chr)] + this_tss
	chunk_manhattan <- ggplot(this_one_chunk_assocs) +
		aes(x=xpos, y=-log10(P), color=factor(chr %% 2)) +
		geom_vline(xintercept=cumul_pos_tss, color="#fd8d3c") +
		geom_point(shape=19) +
		geom_hline(yintercept=-log10(this_prs_assoc_p), linetype=2, color="#7a0177") +
		geom_hline(yintercept=-log10(0.05), linetype=2, color="red") +
		scale_colour_manual(guide=FALSE, values=c("0"="#9ecae1", "1"="#3182bd")) +
		scale_x_continuous(name="Chromosome", breaks=xticks, labels=gsub("21", "", 1:22), expand=c(0,0.0001)) +
		scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0.05, 0.1))) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(),
			panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
			axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
		)

	# Plot change in association P-value as number of chunks is decreased
  this_leaveout_assocs <- leaveout_assocs[PRS == this_prs & Gene == this_gene]
	pval_change <- ggplot(this_leaveout_assocs) +
		aes(x=pct_removed*100, y=-log10(P)) +
		geom_line(colour="#3182bd") +
		geom_hline(yintercept=-log10(0.05), linetype=2, color="red") +
		scale_x_continuous(name="Proportion of genome removed from PRS", limits=c(0, 100), expand=expansion(mult=c(0.01,0.01))) +
		scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0.05, 0.05))) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(),
			panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
			axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
		)
  
  g <- plot_grid(pqtl_manhattan, chunk_manhattan, pval_change, nrow=3)
  ggsave(g, width=7.2, height=3, file=sprintf("analyses/pub/cardiometabolic_proteins/review2/polygenicity_%s_%s.pdf", this_prs, this_gene)) 
}


