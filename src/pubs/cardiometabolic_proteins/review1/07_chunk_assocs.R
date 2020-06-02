library(data.table)
library(foreach)
library(doMC)

registerDoMC(6)

source("src/utilities/prot_pval.R")

out_dir="analyses/pub/cardiometabolic_proteins/review1"

# Load three scores, and filter to set of variants used (ultimately just need the chr and pos)
scores <- rbind(idcol="PRS",
  CAD_metaGRS=fread("data/GRS_resources/CAD_metaGRS/grs_weights.txt.gz"),
  T2D_2018=fread("data/GRS_resources/T2D_2018/grs_weights.txt.gz"),
  CKD_2019=fread("data/GRS_resources/CKD_2019/grs_weights.txt.gz")
)

# Drop strand ambiguous variants:
flip_strand <- function(x) {
  # Swap each letter for a dummy, we need this intermediate
  # step so we can distinguish between alleles when swapping.
  # E.g if we did A -> T then T -> A we'd end up with all A's
  # and no T's. instead we do A -> V -> T and T -> X -> A.
  x <- gsub("A", "V", x)
  x <- gsub("T", "X", x)
  x <- gsub("C", "Y", x)
  x <- gsub("G", "Z", x)
  x <- gsub("V", "T", x)
  x <- gsub("X", "A", x)
  x <- gsub("Y", "G", x)
  x <- gsub("Z", "C", x)
  return(x)
}
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
fwrite(scores, sep="\t", quote=FALSE, file=sprintf("%s/filtered_score_files.txt", out_dir))

# Next we want to load in the GWAS associations for these variants for all 
# proteins of interest
soma_assocs <- fread("analyses/pub/cardiometabolic_proteins/all_assocs.tsv")
soma_assocs <- soma_assocs[Prot.FDR < 0.05]
soma_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

aptamer_set <- soma_assocs[soma_info, on = .(Aptamer=SeqId), nomatch=0, .(SOMAMER_ID, Gene, Target, UniProt)]
aptamer_set <- unique(aptamer_set)

chrpos <- unique(scores[, .(chr, pos)])

gwas <- foreach(gene_id = unique(aptamer_set$Gene), .combine=rbind) %do% { 
  foreach(chr_id = 1:22, .combine=rbind) %do% {
  	all_apt <- foreach(soma_id = unique(aptamer_set[Gene == gene_id, SOMAMER_ID]), .combine=rbind) %do% {
      this_apt <- fread(sprintf("data/full_pQTL_summary_stats/%s/%s_chrom_%s_meta_1.tbl.gz", soma_id, soma_id, chr_id))
      this_apt <- chrpos[chr == chr_id][this_apt, on = .(pos=position), nomatch=0]
      this_apt <- this_apt[, .(chr, pos, effect_allele=toupper(Allele1), other_allele=toupper(Allele2),
                               effect=Effect, se=StdErr, pval=10^`log(P)`)]
      this_apt
    }
    all_apt <- all_apt[, .(effect=mean(effect), se=mean(se), pval=prot_pvalue(pval, effect)), by=.(chr, pos, effect_allele, other_allele)]
    all_apt[, Gene := gene_id]
    return(all_apt)
  }
}

# Add flag for cis or trans regions:
cis_windows <- soma_info[Gene.Name %in% unique(gwas$Gene), .(Gene=Gene.Name, chr=as.integer(chr), TSS=as.integer(start))]
cis_windows[, start := pmax(TSS - 1e6, 0)]
cis_windows[, end := TSS + 1e6]

# Define high complexity regions
# See flashpca exclusion regions: https://github.com/gabraham/flashpca
# coordinates are HG29
complex_ld <- data.table(
  region_chr=c(5, 6, 8, 11),
  region_start=c(44000000, 25000000, 8000000, 45000000),
  region_end=c(51500000, 33500000, 12000000, 57000000),
  region_name=c("r1", "MHC", "r3", "r4")
)


# If 1 MB window inside any region, redefine window as region
cis_windows[complex_ld, on = .(chr = region_chr, start > region_start, end < region_end),
          c("start", "end") := .(region_start, region_end)]
# If 1 MB window overlaps the region, extend the window:
cis_windows[complex_ld, on = .(chr = region_chr, start < region_end, end > region_end),
          c("start") := .(region_start)]
cis_windows[complex_ld, on = .(chr = region_chr, end > region_start, start < region_start),
          c("end") := .(region_end)]

gwas[, cis := FALSE]
gwas[cis_windows, on = .(Gene, chr, pos >= start, pos <= end), cis := TRUE]
fwrite(gwas, sep="\t", quote=FALSE, file=sprintf("%s/sig_prot_score_variant_pQTL_effects.txt", out_dir))
fwrite(cis_windows, sep="\t", quote=FALSE, file=sprintf("%s/cis_regions.txt", out_dir))

# Next, load the protein data and covariates so we can test each 10MB window for protein association.
prot <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")
prot[soma_info, on = .(variable), c("Gene", "Aptamer") := .(Gene.Name, SeqId)]
prot <- prot[Aptamer %in% soma_assocs$Aptamer]
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]
batch <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")
covar <- batch[pcs, on = .(IID), nomatch=0]

# Load the chunked score levels
chunks <- fread("analyses/pub/cardiometabolic_proteins/review1/chunked_scores.sscore.gz")
chunks <- chunks[IID %in% covar$IID]
chunks[, PRS := gsub("_chr.*$", "", score_chunk)]
chunks[, chr := as.integer(gsub("_.*$", "", gsub(".*chr", "", score_chunk)))]
chunks[, chunk_10MB_num := as.integer(gsub(".*_", "", score_chunk))]
chunks[, score_chunk := NULL]
fwrite(chunks, sep="\t", quote=FALSE, file=sprintf("%s/chunk_scores_soma_subset.txt", out_dir))

# Test associations between each chunk and each protein
soma_assocs[PRS == "Chronic Kidney Disease", PRS := "CKD_2019"]
soma_assocs[PRS == "Coronary Artery Disease", PRS := "CAD_metaGRS"]
soma_assocs[PRS == "Type 2 Diabetes", PRS := "T2D_2018"]
pairs <- soma_assocs[, .(PRS, Gene, Aptamer)]

registerDoMC(5)
chunk_assocs <- foreach(prs_id = unique(pairs$PRS), .combine=rbind) %:%
  foreach(chr_id = 1:22, .combine=rbind) %:% 
    foreach(chunk_id = chunks[PRS == prs_id & chr == chr_id, unique(chunk_10MB_num)], .combine=rbind) %:%
			foreach(gene_id = pairs[PRS == prs_id, unique(Gene)], .combine=rbind) %dopar% {
				apt_assocs <- foreach(apt_id = pairs[PRS == prs_id & Gene == gene_id, Aptamer], .combine=rbind) %do% {
          dat <- chunks[PRS == prs_id & chr == chr_id & chunk_10MB_num == chunk_id]
          dat <- dat[prot[Aptamer == apt_id], on = .(IID), nomatch=0]
          dat <- dat[covar, on = .(IID), nomatch=0]
          l1 <- lm(value ~ scale(score_sum) + factor(batch) + PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10, data=dat)
          cf <- coef(summary(l1))
          ci <- confint(l1)
          data.table(Beta=cf[2,1], L95=ci[2,1], U95=ci[2,2], Pvalue=cf[2,4])
        }
        apt_assocs[, .(PRS = prs_id, chr = chr_id, chunk_10MB_num = chunk_id, Gene = gene_id,
                       Beta = mean(Beta), L95=mean(L95), U95=mean(U95), Pvalue=prot_pvalue(Pvalue, Beta))]
}

# Add cis flag
chunk_assocs[, chunk_start := chunk_10MB_num * 10e6]
chunk_assocs[, chunk_end := chunk_start + 10e6 - 1]
chunk_assocs[, cis := FALSE]
chunk_assocs[cis_windows, on = .(Gene, chr, chunk_start >= start, chunk_start <= end), cis := TRUE] # chunk start in window
chunk_assocs[cis_windows, on = .(Gene, chr, chunk_end >= start, chunk_end <= end), cis := TRUE] # chunk end in window
chunk_assocs[cis_windows, on = .(Gene, chr, chunk_start <= start, chunk_end >= end), cis := TRUE] # cis window entirely inside chunk
fwrite(chunk_assocs, sep="\t", quote=FALSE, file=sprintf("%s/one_chunk_assocs.txt", out_dir))

# Now, progressively leave out each chunk from the score (ordered by P-value) to see how
# progressively removing 10 MB chunks affects the scores
chunk_assocs <- chunk_assocs[order(Pvalue)][order(Gene)][order(PRS)]
chunk_assocs[, chunk_rank := seq_len(.N), by=.(PRS, Gene)]

polygenicity <- foreach(prs_id = unique(pairs$PRS), .combine=rbind) %:%
  foreach(gene_id = pairs[PRS == prs_id, unique(Gene)], .combine=rbind) %:%
    foreach(chunk_rank_id = chunk_assocs[PRS == prs_id & Gene == gene_id, chunk_rank], .combine=rbind) %dopar% {
      gc()
      agg_scores <- chunks[chunk_assocs[PRS == prs_id & Gene == gene_id & chunk_rank >= chunk_rank_id], on = .(PRS, chr, chunk_10MB_num)]
      agg_scores <- agg_scores[, .(score_sum=sum(score_sum)), by=.(IID)]
		  apt_assocs <- foreach(apt_id = pairs[PRS == prs_id & Gene == gene_id, Aptamer], .combine=rbind) %do% {
				dat <- agg_scores[prot[Aptamer == apt_id], on = .(IID), nomatch=0]
				dat <- dat[covar, on = .(IID), nomatch=0]
				l1 <- lm(value ~ scale(score_sum) + factor(batch) + PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10, data=dat)
				cf <- coef(summary(l1))
				ci <- confint(l1)
				data.table(Beta=cf[2,1], L95=ci[2,1], U95=ci[2,2], Pvalue=cf[2,4])
			}
			apt_assocs[, .(PRS = prs_id, Gene = gene_id, n_chunks = chunk_assocs[PRS == prs_id & Gene == gene_id, max(chunk_rank)] - chunk_rank_id + 1, 
										 Beta = mean(Beta), L95=mean(L95), U95=mean(U95), Pvalue=prot_pvalue(Pvalue, Beta))]
}
fwrite(polygenicity, sep="\t", quote=FALSE, file=sprintf("%s/progressive_chunk_leaveout_assocs.txt", out_dir))

