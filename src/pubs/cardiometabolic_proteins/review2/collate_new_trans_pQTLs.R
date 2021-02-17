library(data.table)
source("src/utilities/flip_strand.R")

# Load new trans-pQTLs
new_trans <- fread("analyses/pub/cardiometabolic_proteins/review2/new_trans_pQTLs.txt")

# Obtain their MAF
chr_stats <- foreach(chr_idx = 1:22, .combine=rbind) %do% {
	fread(sprintf("data/INTERVAL/reference_files/imputed_genotypes/impute_%s_interval.snpstats", chr_idx))
}
new_trans[chr_stats, on = .(chr=chromosome, pos=position, EA=minor_allele), EAF := MAF]
new_trans[chr_stats, on = .(chr=chromosome, pos=position, EA=major_allele), EAF := 1-MAF]

# Ambiguous SNP MAF cutoff
ambig_max_maf <- 0.42 # set to 0 to drop ambiguous SNPs altogether

# Other options
biallelic_only <- TRUE
snps_only <- TRUE

# Flag ambiguous SNPs for removal
new_trans[, ambig := EA == flip_strand(OA)]
new_trans <- new_trans[!(ambig) | EAF < ambig_max_maf | EAF > 1 - ambig_max_maf]

# Filter to biallelic SNPs if requested
if (biallelic_only) {
  multi <- unique(new_trans[, .(chr, pos, EA, OA)])
  multi <- multi[,.N,by=.(chr, pos)][N > 1]
  new_trans <- new_trans[!multi, on = .(chr, pos)]
}

if (snps_only) {
  new_trans <- new_trans[nchar(EA) == 1 & nchar(OA) == 1]
}

# Load protein information.
prot_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

# Filter to aptamers passing QC
prot_info <- prot_info[Type == "Protein"]

# Select columns
prot_info <- prot_info[, .(SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot,
                           Gene=Gene.Name, chr, start, Cross_Reactivity=Characterization.Info)]

# Fix bad entries (Aptamers for the same target with different/missing gene/uniprot information)
prot_info[Target == "14-3-3 protein family", UniProt := "P61981|Q04917"]
prot_info[Target == "Induced myeloid leukemia cell differentiation protein Mcl-1",
  c("UniProt", "Gene", "chr", "start") :=
  .("Q07820", "MCL1", "1", "150547027")]
prot_info[Target == "Protein delta homolog 1",
  c("UniProt", "Gene", "chr", "start") :=
  .("P80370", "DLK1", "14", "101193202")]
prot_info[Target == "Stromal cell-derived factor 1",
  c("UniProt", "Gene", "chr", "start") :=
  .("P48061", "CXCL12", "10", "44865601")]

# Filter to aptamers with cis-pQTLs:
prot_info <- prot_info[Aptamer %in% unique(new_trans$Aptamer)] # N = 961 aptamers, 890 proteins
new_trans <- new_trans[Aptamer %in% prot_info$Aptamer]

# Filter to lead SNP per 1Mb region
ind_pQTLs <- new_trans[0]
iter = 1
while(nrow(new_trans) > 0) {
  cat(sprintf("Iteration %s, pQTLs remaining: %s\n", iter, new_trans[,.N]))
  top_chr <- new_trans[order(-abs(Effect))][order(P)][order(chr)][order(Aptamer)][, .SD[1], by=.(Aptamer, chr)]
  ind_pQTLs <- rbind(ind_pQTLs, top_chr)

  top_chr <- top_chr[, .(Aptamer, chr, win_start = pos - 1e6, win_end = pos + 1e6)]

  new_trans <- new_trans[!top_chr, on = .(Aptamer, chr, pos >= win_start, pos <= win_end)]
  iter <- iter + 1
}

# Orient effects to minor allele in INTERVAL
ind_pQTLs[EAF > 0.5, c("EA", "OA", "EAF", "Effect") := .(OA, EA, 1-EAF, -Effect)]

# Add protein information
ind_pQTLs <- prot_info[ind_pQTLs, on = .(Target, UniProt, Gene, Aptamer), nomatch=0]

# Select columns
ind_pQTLs <- ind_pQTLs[, .(Target, UniProt, Gene, Aptamer, chr, start, pQTL.chr=i.chr, pQTL.pos=pos, EA,
                           OA, EAF, beta=Effect, SE, P)]

# Order rows
ind_pQTLs <- ind_pQTLs[order(P)][order(pQTL.chr)][order(Aptamer)][order(Target)][order(Gene)]

# Write out
fwrite(ind_pQTLs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/new_trans_pQTLs_curated.txt")

