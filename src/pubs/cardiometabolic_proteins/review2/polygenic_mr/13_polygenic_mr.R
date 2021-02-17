library(data.table)
library(foreach)
library(doMC)
library(openxlsx)
library(MendelianRandomization)
library(ggplot2)
library(ggrastr)

source("src/utilities/flip_strand.R")
source("src/utilities/prot_pval.R")

# Set up parallel environment
if (!exists("ncores")) {
  ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
  ncores <- as.integer(ncores)
  if(is.na(ncores)) ncores <- 1
}
registerDoMC(ncores)
setDTthreads(ncores)

# Load PRS-variants mapped to INTERVAL
scores <- fread("analyses/pub/cardiometabolic_proteins/review2/filtered_score_files.txt")

# Load associations between each PRS 10 Mb chunk and each protein:
one_chunk_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/one_chunk_assocs.txt")

# Load polygenicity metric:
polygenicity <- fread("analyses/pub/cardiometabolic_proteins/review2/polygenicity.txt")

# Map PRS variants to Disease GWAS for MR
prs <- c("CAD_PRS"="CAD_metaGRS", "IS_PRS"="Stroke_metaGRS", "CKD_PRS"="CKD_2019", "T2D_PRS"="T2D_2018")

# Load GWASs of interest
gwass <- sprintf("analyses/mendelian_randomisation/%s/gwas_summary_stats", prs)
gwass <- sapply(gwass, list.files, full.names=TRUE)
if(length(gwass) == 0) stop(sprintf("No gwas summary stats found in analyses/mendelian_randomisation/%s/gwas_summary_stats", prs))
gwass <- data.table(PRS=names(prs), gwas_dir=gwass, gwas_name=basename(gwass))

# Variant filtering options
biallelic_only <- TRUE
snps_only <- TRUE

# Load protein information.
prot_info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

# Filter to aptamers passing QC
prot_info <- prot_info[Type == "Protein"]

# Select columns
prot_info <- prot_info[, .(SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, Gene=Gene.Name)]

# Filter 
prot_info <- prot_info[Target %in% polygenicity$Target]

# Get LD block info
ld_stats <- unique(one_chunk_assocs[, .(ld_block, chr, block_start, block_end, block_size)])
max_chr_block <- ld_stats[,.SD[which.max(ld_block)],by=chr]

# Get GWAS summary statistics mapping to PRS variants
prs_gwas <- foreach(gwas_idx = gwass[,.I], .combine=rbind) %do% {
  this_gwas <- gwass[gwas_idx, gwas_name]
  this_gwas_dir <- gwass[gwas_idx, gwas_dir]
  this_prs <- gwass[gwas_idx, PRS]

  # Get variants in this PRS
  this_score <- scores[PRS == this_prs]

  # Load GWAS summary stats
  gwas_ss <- foreach(chr_idx = 1:22, .combine=rbind) %dopar% {
    chr_ss <- fread(sprintf("%s/chr%s.txt.gz", this_gwas_dir, chr_idx))
    chr_ss <- chr_ss[pos %in% this_score[chr == chr_idx, pos]]
    chr_ss
  }

  # Filter to biallelic SNPs if requested
  if (biallelic_only) {
    multi <- unique(gwas_ss[, .(chr, pos, EA, OA)])
    multi <- multi[,.N,by=.(chr, pos)][N > 1]
    gwas_ss <- gwas_ss[!multi, on = .(chr, pos)]
  }

  if (snps_only) {
    gwas_ss <- gwas_ss[nchar(EA) == 1 & nchar(OA) == 1]
  }

  # (FYI: no strand ambiguous alleles in PRS so don't need to handle here)
  
  # Check for match by effect allele
  gwas_ss[, match := FALSE]
  gwas_ss[this_score, on = .(chr, pos, EA=effect_allele), match := TRUE]
  
  # If effect allele matches non-effect allele in GWAS, swap alleles
  gwas_ss[this_score, on = .(chr, pos, EA=other_allele), 
    c("EA", "OA", "MAF", "effect", "match") :=
    .(OA, EA, 1 - MAF, -effect, TRUE)]

  # If no match, flip allele strand and try again:
  gwas_ss[!(match), c("EA", "OA") := .(flip_strand(EA), flip_strand(OA))]
  gwas_ss[this_score, on = .(chr, pos, EA=effect_allele), match := TRUE]
  gwas_ss[this_score, on = .(chr, pos, EA=other_allele),
    c("EA", "OA", "MAF", "effect", "match") :=
    .(OA, EA, 1 - MAF, -effect, TRUE)]

  # Drop non-matching variants
  gwas_ss <- gwas_ss[(match)]
 
  # Return
  gwas_ss[, .(PRS=this_prs, GWAS=this_gwas, chr, pos, EA, OA, EAF=MAF, effect, se, P)]
} 

# For each PRS to protein association, go through the LD blocks contributing to the
# polygenicity, load the pQTLs and pick the lead SNP in that LD block.
lead_ld_pQTLs <- foreach(test_idx = polygenicity[,.I], .combine=rbind) %dopar% {
  this_prs <- polygenicity[test_idx, PRS]
  this_target <- polygenicity[test_idx, Target]
  this_uniprot <- polygenicity[test_idx, UniProt]
  this_gene <- polygenicity[test_idx, Gene]
  this_aptamers <- prot_info[Target == this_target, SOMAMER_ID]
  this_score <- prs_gwas[PRS == this_prs]

  # Which LD blocks contribute to the polygenicity?
  poly_ld <- one_chunk_assocs[PRS == this_prs & Target == this_target]
  poly_ld <- poly_ld[order(P)]
  poly_ld <- poly_ld[1:polygenicity[test_idx, LD_blocks_removed]]
  poly_ld[max_chr_block, on = .(ld_block), block_end := x.block_end + 1] # if last block on chrosome, add 1 so we can match last pos

  # Load pQTLs for this aptamer in these LD blocks and filter to variants present in PRS
  ld_apt_qtls <- foreach(chr_idx = unique(poly_ld$chr), .combine=rbind) %:% 
    foreach(somaid = this_aptamers, .combine=rbind) %do% {
      apt_qtls <- fread(sprintf("data/full_pQTL_summary_stats/%s/%s_chrom_%s_meta_1.tbl.gz", somaid, somaid, chr_idx))
 
      apt_qtls <- apt_qtls[poly_ld[,.(ld_block, chr, block_start, block_end)], on = .(chromosome=chr, position >= block_start, position < block_end),
                           .(SOMAMER_ID=somaid, chr=chromosome, pos=x.position, ld_block, EA=toupper(Allele1), OA=toupper(Allele2), Effect, SE=StdErr, P=10^(`log(P)`))]
      apt_qtls <- apt_qtls[this_score[,.(chr,pos)], on = .(chr, pos), nomatch=0]
      apt_qtls[this_score, on = .(chr, pos, OA=EA), c("EA", "OA", "Effect") := .(x.OA, x.EA, -Effect)]
      apt_qtls
  }

  # Average across multiple aptamers
  ld_prot_qtls <- ld_apt_qtls[, .(Effect=mean(Effect), SE=mean(SE), P=prot_pvalue(P,Effect)), by=.(chr, pos, ld_block, EA, OA)]

  # Get Best pQTL in PRS per LD block
  best_ld_pqtl <- ld_prot_qtls[,.SD[which.min(P)],by=ld_block]

  # Add back in per-aptamer information
  best_ld_pqtl <- merge(best_ld_pqtl, ld_apt_qtls, by=c("ld_block", "chr", "pos", "EA", "OA"), suffixes=c("", ".aptamer"))

  # Annotate
  info <- data.table(PRS = this_prs, Target = this_target, UniProt = this_uniprot, Gene = this_gene)
  cbind(info, best_ld_pqtl)
}

# Add effect allele frequency information
for (chr_idx in 1:22) {
  snpstats <- fread(sprintf("data/INTERVAL/reference_files/imputed_genotypes/impute_%s_interval.snpstats", chr_idx))
  lead_ld_pQTLs[snpstats, on = .(chr=chromosome, pos=position, EA=minor_allele), EAF := MAF]
  lead_ld_pQTLs[snpstats, on = .(chr=chromosome, pos=position, OA=minor_allele), EAF := 1 - MAF]
}

# Orient to minor allele
lead_ld_pQTLs[EAF > 0.5, c("EA", "OA", "EAF", "Effect", "Effect.aptamer") := .(OA, EA, 1 - EAF, -Effect, -Effect.aptamer)]

# Get aptamer ID
lead_ld_pQTLs[prot_info, on = .(SOMAMER_ID), Aptamer := Aptamer]

# Build table containing the pQTL effects, GWAS effects, and PRS weight
comp <- lead_ld_pQTLs[,.(PRS, Target, UniProt, Gene, chr, pos, ld_block, EA, OA, EAF, pQTL.beta=Effect, pQTL.se=SE, pQTL.pval = P, 
                         Aptamer, aptamer.beta=Effect.aptamer, aptamer.se=SE.aptamer, aptamer.pval = P.aptamer)]
comp[prs_gwas, on = .(PRS, chr, pos, EA), c("GWAS", "GWAS.EAF", "GWAS.logOR", "GWAS.se", "GWAS.pval") := .(GWAS, i.EAF, effect, se, P)]
comp[prs_gwas, on = .(PRS, chr, pos, EA=OA), c("GWAS", "GWAS.EAF", "GWAS.logOR", "GWAS.se", "GWAS.pval") := .(i.GWAS, 1-i.EAF, -effect, se, P)]
comp[,GWAS.pval := as.numeric(GWAS.pval)]
comp[scores, on = .(PRS, chr, pos, EA=effect_allele), PRS.weight := weight]
comp[scores, on = .(PRS, chr, pos, EA=other_allele), PRS.weight := -weight]
comp <- comp[,.(Target, UniProt, Gene, chr, pos, ld_block, EA, OA, EAF, pQTL.beta, pQTL.se, pQTL.pval, GWAS, GWAS.EAF, GWAS.logOR, GWAS.se, GWAS.pval, PRS, PRS.weight, 
                Aptamer, pQTL.beta.aptamer=aptamer.beta, pQTL.se.aptamer=aptamer.se, pQTL.pval.aptamer=aptamer.pval)]
fwrite(comp, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/lead_ldblock_snp.txt")

# Run MR for every test
poly_mr <- foreach(test_idx = polygenicity[,.I], .combine=rbind) %dopar% {
  this_prs <- polygenicity[test_idx, PRS]
  this_target <- polygenicity[test_idx, Target]
  this_uniprot <- polygenicity[test_idx, UniProt]
  this_gene <- polygenicity[test_idx, Gene]

  this_comp <- unique(comp[PRS == this_prs & Target == this_target, .(GWAS, pQTL.beta, pQTL.se, GWAS.logOR, GWAS.se)])
  this_gwas <- this_comp[,unique(GWAS)]

  if (nrow(this_comp) < 3) return(NULL)

  # Construct MR Input object
  mri <- this_comp[, mr_input(bx=pQTL.beta, bxse=pQTL.se, by=GWAS.logOR, byse=GWAS.se)]

  # Run 5 MR methods
  ivw <- mr_ivw(mri)
  simple_median <- mr_median(mri, weighting="simple")
  weighted_median <- mr_median(mri, weighting="weighted")
  weighted_mode <- mr_mbe(mri, weighting="weighted", stderror="simple")
  egger <- mr_egger(mri)

  # Extract relevant coefficients
  extract_coef <- function(mro, egger=FALSE) {
    if (!egger) {
      data.table(mr_estimate = mro@Estimate,
                 mr_se = mro@StdError,
                 mr_L95 = mro@CILower,
                 mr_U95 = mro@CIUpper,
                 mr_pval = mro@Pvalue)
    } else {
      data.table(mr_estimate = c(mro@Estimate, mro@Intercept),
                 mr_se = c(mro@StdError.Int),
                 mr_L95 = c(mro@CILower.Est, mro@CILower.Int),
                 mr_U95 = c(mro@CIUpper.Est, mro@CIUpper.Int),
                 mr_pval = c(mro@Causal.pval, mro@Pleio.pval))
    }
  }

  mro <- rbind(idcol="Method",
    "IVW"=extract_coef(ivw),
    "Median"=extract_coef(simple_median),
    "Weighted Median"=extract_coef(weighted_median),
    "Weighted Mode"=extract_coef(weighted_mode),
    "MR Egger"=extract_coef(egger, TRUE)
  )
  mro[6, Method := "(Intercept)"]

  info <- data.table(PRS = this_prs, GWAS = this_gwas, Target = this_target, UniProt = this_uniprot, Gene = this_gene)
  mro <- cbind(info, mro)
  mro
}

# Take median of causal estimates
mr_summary <- poly_mr[Method != "(Intercept)",
    .(mr_estimate = median(mr_estimate), mr_se = median(mr_se),
      mr_L95 = median(mr_L95), mr_U95 = median(mr_U95),
      mr_pval = prot_pvalue(mr_pval, mr_estimate, median)),
  by = .(PRS, GWAS, Target, UniProt, Gene)]

# Add pleiotropy p-value
mr_summary[poly_mr[Method == "(Intercept)"], on = .(PRS, GWAS, Target, UniProt, Gene), pleiotropy_pval := i.mr_pval]

# Flag pQTL status
sun <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=5, startRow=3)
sun <- as.data.table(sun)
head_row <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, rows=5:6, fillMergedCells=TRUE)
head_row <- gsub(".NA$", "", paste(colnames(head_row), as.vector(head_row[1,]), sep="."))
pQTL_info <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, startRow=6)
colnames(pQTL_info) <- head_row
pQTL_info <- as.data.table(pQTL_info)
sun[pQTL_info, on = .(SOMAmer.ID, Sentinel.variant=`Sentinel.variant*`), type := `cis/.trans`]
sun <- sun[,.(SOMAMER_ID=SOMAmer.ID, rsid=Conditional.variant, type)] 
sun <- sun[SOMAMER_ID %in% prot_info$SOMAMER_ID]
sun <- sun[prot_info, on = .(SOMAMER_ID), nomatch=0]

new_trans <- fread("analyses/pub/cardiometabolic_proteins/review2/new_trans_pQTLs.txt")
new_trans <- new_trans[Target %in% prot_info$Target]
new_trans <- new_trans[,.SD[which.min(P)],by=.(Target, UniProt, Gene)]

cis <- fread("analyses/pub/cardiometabolic_proteins/review2/cis_pQTLs.txt")
cis <- cis[Target %in% prot_info$Target]

mr_summary[, c("trans_pQTLs", "cis_pQTLs", "three_cis") := FALSE]
mr_summary[Target %in% c(sun[type == "trans", Target], new_trans$Target), trans_pQTLs := TRUE]
mr_summary[Target %in% c(sun[type == "cis", Target], cis$Target), cis_pQTLs := TRUE]
mr_summary[cis_mr, on = .(PRS, GWAS, Target), three_cis := TRUE]
mr_summary[!(trans_pQTLs) & !(cis_pQTLs), class := "(0) No pQTLs"]
mr_summary[(trans_pQTLs) & !(cis_pQTLs), class := "(1) Trans pQTLs only"]
mr_summary[(cis_pQTLs), class := "(2) At least one cis-pQTL"]
mr_summary[(three_cis), class := "(3) Three or more cis-pQTLs"]
mr_summary[, c("trans_pQTLs", "cis_pQTLs", "three_cis") := NULL]

# Combine with individual estimates
poly_mr <- merge(mr_summary, poly_mr, by=c("PRS", "GWAS", "Target", "UniProt", "Gene"), suffixes=c(".consensus", ""))

# Order
poly_mr <- poly_mr[order(mr_pval.consensus)][order(pleiotropy_pval < 0.05)][order(PRS)]
fwrite(poly_mr, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/polygenic_mr.txt")

# Add LD-block information to instruments
poly_ivs <- comp
ld_blocks <- fread("data/Berisa_etal_2016/EUR_1000G_ind_ld_blocks.bed")
ld_blocks[, chr := as.integer(gsub("chr", "", chr))]
ld_blocks[, ld_block := .I]
poly_ivs <- poly_ivs[ld_blocks, on = .(ld_block), nomatch=0]

# Write out instruments
poly_ivs <- poly_ivs[unique(poly_mr[,.(PRS, Target, UniProt, Gene)]), on = .(PRS, Target, UniProt, Gene)]
fwrite(poly_ivs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/polygenic_ivs.txt")

################################
# Plot all MR estimates
################################

tests <- unique(poly_mr[,.(PRS, Gene)])
foreach(test_idx = tests[,.I]) %do% {
  this_prs <- tests[test_idx, PRS]
  this_gene <- tests[test_idx, Gene]
 
	this_ivs <- poly_ivs[Gene == this_gene & PRS == this_prs]
	this_ivs[, c("Beta", "Beta.LSE", "Beta.USE", "logOR", "logOR.LSE", "logOR.USE") :=
		.(pQTL.beta, pQTL.beta - pQTL.se, pQTL.beta + pQTL.se,
			GWAS.logOR, GWAS.logOR - GWAS.se, GWAS.logOR + GWAS.se)]

	this_mr <- poly_mr[Gene == this_gene & PRS == this_prs & Method == "(Intercept)",
											 .(logOR=mr_estimate.consensus, L95=mr_L95.consensus, U95=mr_U95.consensus,
												 int=mr_estimate, int.L95=mr_L95, int.U95=mr_U95)]

	rdt <- this_mr[, .(x = c(-1, 0, 1),
		ymin = c(int.L95 - U95, int.L95, int.L95 + L95),
		ymax = c(int.U95 - L95, int.U95, int.U95 + U95))]
	plotlim <- this_ivs[, .(x=c(-1, 0, 1), xmult=c(abs(min(Beta.LSE)), 0, max(Beta.USE)))]
	expand <- plotlim[, (max(xmult) - min(xmult))*0.05]
	xmult <- plotlim[,xmult + expand]
	ymult <- c(xmult[1], 1, xmult[3])
	rdt[, c("x", "ymin", "ymax") := .(x*xmult, ymin * ymult, ymax * ymult)]

	rdt[, ymin := exp(ymin)]
	rdt[, ymax := exp(ymax)]

	g <- ggplot(this_ivs) +
		aes(x = Beta, xmin = Beta.LSE, xmax = Beta.USE,
				y = exp(logOR), ymin = exp(logOR.LSE), ymax = exp(logOR.USE)) +
		geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#ffeda0") +
		geom_abline(intercept=exp(this_mr$int), slope=exp(this_mr$logOR) - 1, linetype=2, color="#f03b20") +
		geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
		geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
		geom_errorbarh(height=0, size=0.5, color="black", alpha=0.7) +
		geom_errorbar(width=0, size=0.5, color="black", alpha=0.7) +
		geom_point(shape = 19, size=1, color="black") +
		scale_x_continuous(name = "SD effect on Protein", expand=c(0,0)) +
		scale_y_continuous(name = "Odds Ratio") +
		theme_bw() + theme(
			axis.title=element_text(size=10), axis.text=element_text(size=8),
			panel.grid=element_blank(), legend.position="bottom"
		)
	ggsave(g, width=3.4, height=3.4, useDingbats=FALSE,
         file=sprintf("analyses/pub/cardiometabolic_proteins/review2/%s_%s_polygenic_mr.pdf", this_prs, this_gene))
}

#################################################
# Compare estimated effects to PRS associations
#################################################

prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])
prs_assocs <- prs_assocs[PRS.FDR < 0.05]

prs_assocs[polygenicity, on = .(PRS, Gene), polygenicity := round(pct_removed*100)]
prs_assocs <- prs_assocs[polygenicity > 0]

consensus_mr <- poly_mr[Method == "(Intercept)", .(PRS, Gene, logOR=mr_estimate.consensus, 
                        L95=mr_L95.consensus, U95=mr_U95.consensus, P=mr_pval.consensus,
                        int=mr_estimate, int.L95=mr_L95, int.U95=mr_U95, int.pval=mr_pval)]

prs_assocs <- prs_assocs[consensus_mr, on = .(PRS, Gene), nomatch=0]
prs_assocs[, sig_mr := ifelse(P < 0.05 & int.pval > 0.05, TRUE, FALSE)]

g <- ggplot(prs_assocs) +
  aes(x = PRS.Beta, xmin = PRS.L95, xmax = PRS.U95,
      y = exp(logOR), ymin = exp(L95), ymax = exp(U95),
      color = sig_mr) +
  geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.8, size=0.5) +
  geom_errorbar(width=0, alpha=0.8, size=0.5) +
  geom_point(shape = 19, size=1.3) +
  facet_wrap(~ PRS, nrow=1) +
  scale_x_continuous(name = "PRS effect (95% CI)") +
  scale_y_continuous(name = "Causal OR (95% CI)") +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g, width=7, height=2.6, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/PRS_polygenic_MR_compare.pdf")

# And to cis-PQTL causal estimates
cis_mr <- fread("analyses/pub/cardiometabolic_proteins/review2/collated_mr.txt")
cis_mr <- unique(cis_mr[,.(GWAS, Target, UniProt, Gene, MR=exp(mr_estimate), MR.L95=exp(mr_L95), MR.U95=exp(mr_U95), MR.P=mr_pval,
                                 MR.FDR=mr_fdr, MR.pleiotropy=pleiotropy_pval, colocalization)])
cis_mr[GWAS == "CAD", PRS := "CAD_PRS"]
cis_mr[GWAS == "CKD", PRS := "CKD_PRS"]
cis_mr[GWAS == "StrokeIS", PRS := "IS_PRS"]
cis_mr[GWAS == "T2DadjBMI", PRS := "T2D_PRS"]
prs_assocs <- prs_assocs[cis_mr, on = .(PRS, Target, UniProt, Gene), nomatch=0]
prs_assocs[, sig_cis := ifelse(MR.P < 0.05 & MR.pleiotropy > 0.05, TRUE, FALSE)]
prs_assocs[, anno := "no"]
prs_assocs[(sig_mr), anno := "poly"]
prs_assocs[(sig_cis), anno := "cis"]
prs_assocs[(sig_cis) & (sig_mr), anno := "both"]

g <- ggplot(prs_assocs) +
  aes(x = MR, xmin = MR.L95, xmax = MR.U95,
      y = exp(logOR), ymin = exp(L95), ymax = exp(U95),
      color = anno) +
  geom_hline(yintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=1, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.8, size=0.5) +
  geom_errorbar(width=0, alpha=0.8, size=0.5) +
  geom_point(shape = 19, size=1.3) +
  facet_wrap(~ PRS, nrow=1) +
  scale_x_continuous(name = "Causal OR (95% CI)") +
  scale_y_continuous(name = "Causal OR (95% CI)") +
  theme_bw() + theme(
    axis.title=element_text(size=8), axis.text=element_text(size=8),
    panel.grid=element_blank(), legend.position="bottom"
  )
ggsave(g, width=7, height=2.6, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/MR_compare.pdf")




###########################################
# Plot CRYZL1 on CAD instrument selection
###########################################

leaveout_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/leaveout_chunk_assocs.txt")
leaveout_assocs <- leaveout_assocs[PRS == "CAD_PRS" & Gene == "CRYZL1"]

p1 <- ggplot(leaveout_assocs) +
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
ggsave(p1, width=3.6, height=1.1, file="analyses/pub/cardiometabolic_proteins/review2/crzyl1_polygenicity.pdf")

one_chunk_assocs <- one_chunk_assocs[PRS == "CAD_PRS" & Gene == "CRYZL1"][order(P)]
one_chunk_assocs[, poly := FALSE]
one_chunk_assocs[1:polygenicity[Gene == "CRYZL1", LD_blocks_removed], poly := TRUE]
one_chunk_assocs[, block_idx := .I]

p2 <- ggplot(one_chunk_assocs[P < 0.05]) +
    aes(x=block_idx, y=-log10(P), color=poly) +
    geom_point(shape=19, size=0.8) +
    scale_colour_manual(name="", guide=FALSE, values=c("FALSE"="#636363", "TRUE"="#3182bd")) +
    scale_x_continuous(name="LD blocks in order of removal from CAD PRS", expand=expansion(mult=c(0.01,0.01))) +
    scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0.05, 0.05))) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
      axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
    )
ggsave(p2, width=3.6, height=1.1, file="analyses/pub/cardiometabolic_proteins/review2/crzyl1_leaveout_blocks.pdf")

cryzl1_pqtls <- fread("data/full_pQTL_summary_stats/CRYZL1.9207.60.3/CRYZL1.9207.60.3_chrom_8_meta_1.tbl.gz")
cryzl1_pqtls <- cryzl1_pqtls[position >= 19492840 & position < 20060856]
cryzl1_pqtls <- cryzl1_pqtls[,.(chr=chromosome, pos=position, A1=toupper(Allele1), A2=toupper(Allele2), P=10^`log(P)`)]
cryzl1_pqtls[, type := "all"]

snpstats <- fread("data/INTERVAL/reference_files/imputed_genotypes/impute_8_interval.snpstats")
cryzl1_pqtls[snpstats, on = .(pos=position), MAF := i.MAF]
cryzl1_pqtls[, gwas_match := FALSE]
cryzl1_pqtls[, drop := FALSE]
cryzl1_multi <- cryzl1_pqtls[,.N,by=pos][N > 1]
cryzl1_pqtls[cryzl1_multi, on = .(pos), drop := TRUE]
cryzl1_pqtls[nchar(A1) > 1 | nchar(A2) > 1, drop := TRUE]
cryzl1_pqtls[!(drop) & flip_strand(A1) == A2 & MAF > 0.42 & MAF < 1 - 0.42, drop := TRUE]

gwas_snps <- fread("analyses/mendelian_randomisation/CAD_metaGRS/gwas_summary_stats/CAD/chr8.txt.gz")
gwas_multi <- gwas_snps[,.N,by=pos][N > 1]
gwas_snps <- gwas_snps[!gwas_multi, on = .(pos)]
gwas_snps <- gwas_snps[nchar(EA) == 1 & nchar(OA) == 1]
gwas_snps[, ambig := EA == flip_strand(OA)]
gwas_snps <- gwas_snps[!(ambig) | MAF < 0.42 | MAF > 1 - 0.42]

cryzl1_pqtls[gwas_snps, on = .(pos, A1=EA), gwas_match := TRUE]
cryzl1_pqtls[gwas_snps, on = .(pos, A1=OA), gwas_match := TRUE]
gwas_snps[, c("EA", "OA") := .(flip_strand(EA), flip_strand(OA))]
cryzl1_pqtls[gwas_snps, on = .(pos, A1=EA), gwas_match := TRUE]
cryzl1_pqtls[gwas_snps, on = .(pos, A1=OA), gwas_match := TRUE]
cryzl1_pqtls[(gwas_match) & !(drop), type := "GWAS"]

cad_prs <- fread("data/GRS_resources/CAD_metaGRS/grs_weights.txt.gz")
cad_prs <- cad_prs[chr == 8]
cryzl1_pqtls[!(drop) & flip_strand(A1) == A2, drop := TRUE]
cryzl1_pqtls[, in_prs := FALSE]
cryzl1_pqtls[cad_prs, on = .(pos, A1=effect_allele), in_prs := TRUE]
cryzl1_pqtls[cad_prs, on = .(pos, A1=other_allele), in_prs := TRUE]
cad_prs[, c("effect_allele", "other_allele") := .(flip_strand(effect_allele), flip_strand(other_allele))]
cryzl1_pqtls[cad_prs, on = .(pos, A1=effect_allele), in_prs := TRUE]
cryzl1_pqtls[cad_prs, on = .(pos, A1=other_allele), in_prs := TRUE]
cryzl1_pqtls[type == "GWAS" & !(drop) & (in_prs), type := "PRS"]

lead <- cryzl1_pqtls[type == "PRS"][order(P)][1]
cryzl1_pqtls[lead, on = .(pos), type := "IV"]

p3 <- ggplot(cryzl1_pqtls) +
    aes(x=pos, y=-log10(P)) +
    geom_point_rast(data=cryzl1_pqtls[type == "all"], shape=19, size=0.5, color="#333333", raster.dpi=600, raster.width=4.8, raster.height=1.1) +
    geom_point_rast(data=cryzl1_pqtls[type == "GWAS"], shape=19, size=0.5, color="#08519c", raster.dpi=600, raster.width=4.8, raster.height=1.1) +
    geom_point_rast(data=cryzl1_pqtls[type == "PRS"], shape=19, size=0.5, color="#542788", raster.dpi=600, raster.width=4.8, raster.height=1.1) +
    geom_point(data=cryzl1_pqtls[type == "IV"], shape=19, size=1, color="#e08214") +
    scale_x_continuous(name="Chromosome 8", expand=expansion(mult=c(0.01,0.01))) +
    scale_y_continuous(name="-log10 P-value", expand=expansion(mult=c(0.05, 0.15))) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
      axis.title = element_text(size=8), axis.text=element_text(size=8), title=element_text(size=10)
    )
ggsave(p3, width=3.6, height=1.1, file="analyses/pub/cardiometabolic_proteins/review2/crzyl1_LD_block_lead_snp.pdf")

