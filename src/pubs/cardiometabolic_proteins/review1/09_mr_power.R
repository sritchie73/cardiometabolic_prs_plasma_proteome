library(data.table)
library(foreach)
library(ggplot2)
library(doMC)
library(openxlsx)
library(cowplot)
registerDoMC(32)
source("src/07_job_scripts/07_helpers/mr_functions.R")
source("src/utilities/prot_pval.R")
source("src/utilities/flip_strand.R")
source("src/utilities/format_pval.R")

out_dir <- "analyses/pub/cardiometabolic_proteins/review1"

instruments <- rbind(
  fread("analyses/mendelian_randomisation/CAD_metaGRS/instruments.tsv"),
  fread("analyses/mendelian_randomisation/T2D_2018/instruments.tsv"),
  fread("analyses/mendelian_randomisation/CKD_2019/instruments.tsv")
)

# Drop variants that were not in the GWAS
instruments <- instruments[!is.na(P.gwas)]

# Drop ApoE - any causal effect likely reflects its
# cross-reactivity to the ApoE2 and ApoE4 isoform aptamers,
# which we know will have a causal effect due to their
# effects on LDL levels
instruments <- instruments[Gene != "APOE"]

# Require proteins to have at least three pQTLs, at least one cis
iv_stats <- instruments[, .(.N, n_cis=sum(type=="cis")), by=.(PRS, Gene)]
instruments <- instruments[iv_stats[N >= 3 & n_cis >= 1], on = .(PRS, Gene)]

# Load protein levels, annotate, and filter
prot <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")
info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
prot[info, on = .(variable), Gene := Gene.Name]
prot <- prot[Gene %in% unique(instruments$Gene)]

# Extract dosages of instruments:
fwrite(unique(prot[,.(IID, IID)]), col.names=FALSE, sep=" ", file=sprintf("%s/somalogic_samples.txt", out_dir))
dmerge <- function(x, y) { merge(x, y, by="IID") }
dosages <- foreach(chr_idx = instruments[,unique(chr)], .combine=dmerge) %do% {
  fwrite(instruments[chr == chr_idx, .(RSID)], col.names=FALSE, file=sprintf("%s/instruments_chr%s.txt", out_dir, chr_idx))

  cmd <- "plink2"
  cmd[2] <- sprintf("--pfile data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/plink_format/pgen/impute_dedup_%s_interval", chr_idx)
  cmd[3] <- sprintf("--extract %s/instruments_chr%s.txt", out_dir, chr_idx)
  cmd[4] <- sprintf("--keep %s/somalogic_samples.txt", out_dir)
  cmd[5] <- "--export A"
  cmd[6] <- sprintf("--out %s/instruments_chr%s_dosages.txt", out_dir, chr_idx)
  cmd[7] <- sprintf("--threads 1 --memory 8192")
  system(paste(cmd, collapse=" "), wait=TRUE)
  
  dosages <- fread(sprintf("%s/instruments_chr%s_dosages.txt.raw", out_dir, chr_idx))
  dosages[, c("FID", "PAT", "MAT", "SEX", "PHENOTYPE") := NULL]
  setnames(dosages, gsub("_.*", "", names(dosages)))
  dosages
}

# Fit multivariable model for each aptamer, and obtain the r2's
r2 <- foreach(grs = unique(instruments$PRS), .combine=rbind) %do% {
  foreach(aptamer = prot[Gene %in% instruments[PRS == grs, unique(Gene)], unique(variable)], .combine=rbind) %do% {
    gene <- info[variable == aptamer, Gene.Name]

    dat <- prot[variable == aptamer]
    dat <- dat[dosages, on = .(IID), nomatch=0]

    model <- sprintf("value ~ %s", paste(instruments[PRS == grs & Gene == gene, RSID], collapse=" + "))

    l1 <- lm(as.formula(model), data=dat)

    data.table(PRS = grs, Aptamer = info[variable == aptamer, SeqId], Gene = gene, r2=summary(l1)$r.squared, adj.r2=summary(l1)$adj.r.squared)
  }
}

# Average estimates for aptamers for SHBG and WFIKKN2
r2 <- r2[,.(r2=mean(r2), adj.r2=mean(adj.r2)), by=.(PRS, Gene)]
r2 <- r2[order(-r2)]

# Investigate overlap of cis-pQTLs with trans-signals
apts <- list.files(path="data/full_pQTL_summary_stats/", full.names=TRUE)
all_pqtl <- foreach(chr_idx = instruments[type == "cis", unique(chr)], .combine=rbind) %do% {
  foreach(apt = apts, .combine=rbind) %dopar% {
    dt = fread(sprintf("%s/%s_chrom_%s_meta_1.tbl.gz", apt, basename(apt), chr_idx))
    dt <- dt[unique(instruments[type == "cis" & chr == chr_idx, .(RSID, chromosome=chr, position=pos)]), on = .(chromosome, position)]
    dt <- dt[,.(Aptamer=basename(apt), RSID, chromosome, position, pval=10^`log(P)`)]
    dt[pval < 0.05] 
  }
}
all_pqtl[info, on = .(Aptamer=SOMAMER_ID), c("Target", "Gene") := .(TargetFullName, Gene.Name)] # Add unique targets
all_pqtl <- all_pqtl[!instruments[type == "cis", .(RSID, Gene)], on = .(RSID, Gene)] # Filter to trans-signals

# Add list of trans-signals to instruments
instruments[all_pqtl[,.N, by=RSID], on = .(RSID), N_trans_0.05 := i.N]
instruments[is.na(N_trans_0.05), N_trans_0.05 := 0]

instruments[all_pqtl[pval < 5e-8,.N, by=RSID], on = .(RSID), N_trans_5e8 := i.N]
instruments[is.na(N_trans_5e8), N_trans_5e8 := 0]

instruments[all_pqtl[pval < 1.5e-11,.N, by=RSID], on = .(RSID), N_trans_gwas := i.N]
instruments[is.na(N_trans_gwas), N_trans_gwas := 0]

plieotropy <- unique(instruments[,.(RSID, Gene, N_trans_0.05, N_trans_5e8, N_trans_gwas)])
plieotropy <- melt(plieotropy, id.vars=c("RSID", "Gene"), variable.name = "threshold", value.name = "n_trans")
plieotropy[, threshold := gsub("N_trans_", "", threshold)]

# rerun MR excluding plieotropic variants
instruments2 <- instruments[N_trans_5e8 == 0 & type == "cis"]
instruments2 <- instruments2[instruments2[,.N,by=.(PRS, Gene)][N >= 3], on = .(PRS, Gene)]

mr2 <- foreach(grs = unique(instruments2$PRS), .combine=rbind) %do% {
  foreach(aptamer = instruments2[PRS == grs, unique(SOMAMER_ID)], .combine=rbind) %do% {

    sub_inst <- instruments2[PRS == grs & SOMAMER_ID == aptamer]

    # Construct MR input object:
    mri <- mr_input(
      bx = sub_inst$Prot.Effect.pQTL,
      bxse = sub_inst$Prot.SE.pQTL,
      by = sub_inst$effect.gwas,
      byse = sub_inst$se.gwas,
      snps = sub_inst$ALT_ID,
      effect_allele = sub_inst$EA,
      other_allele = sub_inst$OA,
      eaf = sub_inst$EAF.pqtl,
      exposure = unique(sub_inst$Gene),
      outcome = grs)

    # Run all available MR methods:
    mr_results <- mr_tryall(mri)

    mr_results[, PRS := grs]
    mr_results[, Aptamer := aptamer]
    mr_results[, Gene := unique(sub_inst$Gene)]
   
    mr_results
  }
}
mr2 <- rbind(
  mr2[method %in% c("IVW", "Simple median", "Weighted median", "Weighted mode (simple SE)")],
  mr2[c(which(method == "MR-Egger"), which(method == "MR-Egger") + 1)]
)
mr2[method == "Weighted mode (simple SE)", method := "Weighted mode"]	
mr2 <- mr2[order(Gene)][order(PRS)]

mr2 <- mr2[,.(mr_estimate = mean(mr_estimate), mr_se = mean(mr_se), mr_L95 = mean(mr_L95), mr_U95 = mean(mr_U95),
              mr_pval = prot_pvalue(mr_pval, mr_estimate)), by=.(PRS, Gene, method)]

mr2_summary <- mr2[method != "(intercept)",
                 .(med_estimate=median(mr_estimate),
                   med_L95=median(mr_L95),
                   med_U95=median(mr_U95),
                   med_pval=median(mr_pval)),
                by = .(PRS, Gene)]
mr2_summary[mr2[method == "(intercept)"], on = .(PRS, Gene), horiz_pleio := mr_pval]

# How many proteins overall can be instrumented?
info <- info[Type == "Protein", .(variable, SeqId, SOMAMER_ID, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, 
                                  Gene=Gene.Name, chr=chr, TSS=start)]
cis <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
cis <- cis[flip_strand(EA) != flip_strand(OA)]

load_gwas <- function(path) {
  foreach(chr = 1:22, .combine=rbind) %dopar% {
    fread(sprintf("%s/chr%s.txt.gz", path, chr))
  }
}

gwas <- rbind(idcol="gwas",
  CAD = load_gwas("analyses/mendelian_randomisation/CAD_metaGRS/gwas_summary_stats/CAD/"),
  T2D = load_gwas("analyses/mendelian_randomisation/T2D_2018/gwas_summary_stats/T2DadjBMI/"),
  CKD = load_gwas("analyses/mendelian_randomisation/CKD_2019/gwas_summary_stats/CKD")
)

# One protein is missing from pre-calculated LD, PRSS2, so compute LD
prot_vars <- cis[SOMAMER_ID == "PRSS2.5034.79.1", .(chr, pos)]
chr_stats <- fread("data/INTERVAL/reference_files/imputed_genotypes/impute_7_interval.snpstats")
prot_vars <- chr_stats[prot_vars, on = .(chromosome = chr, position=pos)]
prot_vars <- prot_vars[,.(SNPID=RSID, RSID=RSID, chromosome, position, A_allele, B_allele)]
prot_vars[, chromosome := paste0("0", chromosome)]
fwrite(prot_vars, sep=" ", quote=FALSE, file=sprintf("%s/PRSS2.5034.79.1_7.vars.txt", out_dir))

aptamer <- "PRSS2.5034.79.1"
chr <- 7
window_start <- info[SOMAMER_ID == aptamer, as.integer(TSS) - 1e6]
window_end <- info[SOMAMER_ID == aptamer, as.integer(TSS) + 1e6]

# First use bgenix to extract 1MB window, to reduce memory usage and
# speed up variant extract by qctool2
cmd <- "bgenix"
cmd <- paste(cmd, sprintf("-g data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/impute_%s_interval.bgen", chr))
if (chr < 10) { # %01d doesnt seem to work here??
	cmd <- paste(cmd, sprintf("-incl-range 0%s:%s-%s", chr, window_start, window_end))
} else {
	cmd <- paste(cmd, sprintf("-incl-range %s:%s-%s", chr, window_start, window_end))
}
cmd <- paste(cmd, sprintf("> %s/%s_%s.cis_window.bgen", out_dir, aptamer, chr))
system(command=cmd, wait=TRUE)

# Run qctool2 command extracting the cis-pQTL variants into new BGEN file
cmd <- "qctool2"
cmd <- paste(cmd, sprintf("-g %s/%s_%s.cis_window.bgen", out_dir, aptamer, chr))
cmd <- paste(cmd, sprintf("-og %s/%s_%s.bgen", out_dir, aptamer, chr))
cmd <- paste(cmd, sprintf("-incl-variants %s/%s_%s.vars.txt", out_dir, aptamer, chr))
system(cmd, wait=TRUE)

# Then use LDstore to calculate LD between variants
cmd <- "ldstore"
cmd <- paste(cmd, sprintf("--bcor %s/%s_%s.bcor", out_dir, aptamer, chr))
cmd <- paste(cmd, sprintf("--bgen %s/%s_%s.bgen", out_dir, aptamer, chr))
cmd <- paste(cmd, "--n-threads 1")
system(command=cmd, wait=TRUE)

# Convert LDstore bcor binary into usable table
cmd <- "ldstore"
cmd <- paste(cmd, sprintf("--bcor %s/%s_%s.bcor_1", out_dir, aptamer, chr))
cmd <- paste(cmd, sprintf("--table %s/%s_%s.ld", out_dir, aptamer, chr))
system(command=cmd, wait=TRUE)

# Clean up
system(sprintf("rm %s/%s_%s.bgen", out_dir, aptamer, chr), wait=TRUE)
system(sprintf("rm %s/%s_%s.cis_window.bgen", out_dir, aptamer, chr), wait=TRUE)
system(sprintf("rm %s/%s_%s.vars.txt", out_dir, aptamer, chr), wait=TRUE)
system(sprintf("rm %s/%s_%s.bcor_1", out_dir, aptamer, chr), wait=TRUE)

# Move LD file to designated directory
system(sprintf("mv %s/%s_%s.ld %s/%s.ld", 
  out_dir, aptamer, chr,
  "analyses/mendelian_randomisation/pqtls/SOMAMERs_cis_LD", aptamer), 
  wait=TRUE)

# Which pQTLs have gwas snps
cis[gwas[gwas == "CAD"], on = .(chr, pos, EA=EA), CAD_GWAS := TRUE]
cis[gwas[gwas == "CAD"], on = .(chr, pos, EA=OA), CAD_GWAS := TRUE]
cis[is.na(CAD_GWAS), CAD_GWAS := FALSE]

cis[gwas[gwas == "T2D"], on = .(chr, pos, EA=EA), T2D_GWAS := TRUE]
cis[gwas[gwas == "T2D"], on = .(chr, pos, EA=OA), T2D_GWAS := TRUE]
cis[is.na(T2D_GWAS), T2D_GWAS := FALSE]

cis[gwas[gwas == "CKD"], on = .(chr, pos, EA=EA), CKD_GWAS := TRUE]
cis[gwas[gwas == "CKD"], on = .(chr, pos, EA=OA), CKD_GWAS := TRUE]
cis[is.na(CKD_GWAS), CKD_GWAS := FALSE]

cis <- cis[CAD_GWAS | T2D_GWAS | CKD_GWAS]

# Determine the independent cis-pQTLs we would use for each MR test
ind_ivs <- function(cis) {
  foreach(apt = unique(cis$SOMAMER_ID), .combine=rbind) %dopar% {
    gc()
    apt_cis <- cis[SOMAMER_ID == apt][order(P)]
    if (nrow(apt_cis) == 1) {
      return(apt_cis[, .(SOMAMER_ID, chr, pos)])
    }
  
    apt_ld <- fread(sprintf("analyses/mendelian_randomisation/pqtls/SOMAMERs_cis_LD/%s.ld", apt))
    apt_ind_ivs <- NULL

    while(nrow(apt_cis) > 0) {
      top <- apt_cis[1]
      apt_cis <- apt_cis[-1]
      inld <- unique(rbind(
        apt_ld[chromosome == top$chr & position1 == top$pos & correlation^2 >= 0.1, .(chr=chromosome, pos=position2)],
        apt_ld[chromosome == top$chr & position2 == top$pos & correlation^2 >= 0.1, .(chr=chromosome, pos=position1)]
      ))
      apt_cis <- apt_cis[!inld, on = .(chr, pos)]
      apt_ind_ivs <- rbind(apt_ind_ivs, top)
    }
    
    apt_ind_ivs[, .(SOMAMER_ID, chr, pos)]
  }
}

rm(gwas, chr_stats); gc()
cis_ind <- rbind(idcol="gwas", 
  CAD=ind_ivs(cis[(CAD_GWAS)]),
  T2D=ind_ivs(cis[(T2D_GWAS)]),
  CKD=ind_ivs(cis[(CKD_GWAS)])
)

n_cis_per_aptamer <- cis_ind[,.N, by=.(gwas, SOMAMER_ID)]
n_cis_per_aptamer <- n_cis_per_aptamer[info, on = .(SOMAMER_ID), nomatch=0]
sufficient_ivs <- unique(n_cis_per_aptamer[N >= 3, .(gwas, Target, UniProt, Gene)])

# Compute r2 for all allelic series
dosages <- foreach(chr_idx = 1:22, .combine=dmerge) %do% {
  pvar <- fread(sprintf("data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/plink_format/pgen/impute_dedup_%s_interval.pvar", chr_idx))
  cis_ind[pvar, on = .(chr=`#CHROM`, pos=POS), rsid := i.ID]
  fwrite(unique(cis_ind[chr == chr_idx, .(rsid)]), col.names=FALSE, file=sprintf("%s/instruments_chr%s.txt", out_dir, chr_idx))

  cmd <- "plink2"
  cmd[2] <- sprintf("--pfile data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/plink_format/pgen/impute_dedup_%s_interval", chr_idx)
  cmd[3] <- sprintf("--extract %s/instruments_chr%s.txt", out_dir, chr_idx)
  cmd[4] <- sprintf("--keep %s/somalogic_samples.txt", out_dir)
  cmd[5] <- "--export A"
  cmd[6] <- sprintf("--out %s/instruments_chr%s_dosages.txt", out_dir, chr_idx)
  cmd[7] <- sprintf("--threads 1 --memory 8192")
  system(paste(cmd, collapse=" "), wait=TRUE)

  dosages <- fread(sprintf("%s/instruments_chr%s_dosages.txt.raw", out_dir, chr_idx))
  dosages[, c("FID", "PAT", "MAT", "SEX", "PHENOTYPE") := NULL]
  dosages
}

# Fit multivariable model for each aptamer, and obtain the r2's
prot <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")
all_r2 <- foreach(gwas_id = c("CAD", "CKD", "T2D"), .combine=rbind) %do% {
  foreach(aptamer = cis_ind[gwas == gwas_id, unique(SOMAMER_ID)], .combine=rbind, .errorhandling="remove") %dopar% {
    apt_target <- info[SOMAMER_ID == aptamer, .(variable, SOMAMER_ID, Target, UniProt, Gene)]

    dat <- prot[variable == apt_target$variable]
    dat <- dat[dosages[,na.omit(c(1, pmatch(cis_ind[gwas == gwas_id & SOMAMER_ID == aptamer, paste0(rsid, "_")], colnames(dosages), duplicates.ok=TRUE))), with=FALSE], on = .(IID)]
    setnames(dat, gsub(":", "_", names(dat)))
    setnames(dat, gsub("^([0-9])", "rs\\1", names(dat)))

    model <- paste("value ~", paste(names(dat)[4:ncol(dat)], collapse=" + "))

    l1 <- lm(as.formula(model), data=dat)

    apt_target[, gwas := gwas_id]
    apt_target[, r2 := summary(l1)$r.squared]
    apt_target[, adj.r2 := summary(l1)$adj.r.squared]

    apt_target
  }
}

# Average estimates for proteins across aptamers 
all_r2_by_gwas_gene <- all_r2[,.(r2=mean(r2), adj.r2=mean(adj.r2)), by=.(gwas, Target, UniProt, Gene)]
all_r2_by_gwas_gene <- all_r2_by_gwas_gene[order(-r2)]

all_r2_by_gene <- all_r2[,.(r2=mean(r2), adj.r2=mean(adj.r2)), by=.(Target, UniProt, Gene)]
all_r2_by_gene <- all_r2_by_gene[order(-r2)]

g <- ggplot(all_r2, aes(x=r2)) +
  geom_density() +
  xlab("r2 explained by cis-pQTLs ld-thinned at r2 > 0.8") +
  theme_bw()
ggsave(g, file=sprintf("%s/all_proteins_cis_pqtl_allelic_series_r2.png", out_dir)) 

sufficient_r2 <- all_r2_by_gene[unique(sufficient_ivs[, .(Target, UniProt, Gene)]), on = .(Target, UniProt, Gene)]
g <- ggplot(sufficient_r2, aes(x=r2)) +
  geom_density() +
  xlab("r2 explained by cis-pQTLs ld-thinned at r2 > 0.8") +
  theme_bw()
ggsave(g, file=sprintf("%s/all_proteins_sufficient_instruments_cis_pqtl_allelic_series_r2.png", out_dir))

cutoffs <- foreach(cut = c(0.01, 0.05, 0.1, 0.25, 0.5), .combine=rbind) %do% {
  sufficient_r2[, .(threshold=sprintf("Variance explained (r2) < %s%%", cut*100), N = sum(r2 < cut), pct = paste0(round(sum(r2 < cut)/.N*100), "%"))]
}

# save all outputs
save(cis_ind, sufficient_ivs, all_r2, sufficient_r2,  
     all_pqtl, mr2, r2, instruments, plieotropy, 
     file=sprintf("%s/mr_power.rda", out_dir))

