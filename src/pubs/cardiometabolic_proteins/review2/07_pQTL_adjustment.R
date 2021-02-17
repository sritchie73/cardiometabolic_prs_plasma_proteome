library(data.table)
library(foreach)
library(openxlsx)
library(doMC)
library(ggplot2)
library(ggrastr)

source("src/utilities/prot_pval.R")

# Set up parallel environment
if (!exists("ncores")) {
  ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
  ncores <- as.integer(ncores)
  if(is.na(ncores)) ncores <- 1
}
registerDoMC(ncores)
setDTthreads(ncores)

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
                                    Gene=Gene.Name, chr, start, INTERVAL_Batch1_SOMALOGIC_QC, INTERVAL_Batch2_SOMALOGIC_QC, INTERVAL_Batch1_CV,
                                    INTERVAL_Batch2_CV)]

# Flag aptamers whose pQTLs have not been mapped by Sun et al.
sinfo[, trans_mapped := INTERVAL_Batch1_SOMALOGIC_QC == "PASS" & INTERVAL_Batch2_SOMALOGIC_QC == "PASS" &
                        INTERVAL_Batch1_CV <= 20 & INTERVAL_Batch2_CV <= 20]

# Load PRS to protein associations
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- prs_assocs[FDR < 0.05]

#################################################################
# Load  pQTLs reported in Sun et al. 2018 for aptamers of interest
#################################################################

# Load published pQTL information
pQTLs <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=5, startRow=3)
pQTLs <- as.data.table(pQTLs)

# Add information about whether each variant is cis or trans
head_row <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, rows=5:6, fillMergedCells=TRUE)
head_row <- gsub(".NA$", "", paste(colnames(head_row), as.vector(head_row[1,]), sep="."))
pQTL_info <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, startRow=6)
colnames(pQTL_info) <- head_row
pQTL_info <- as.data.table(pQTL_info)
pQTLs[pQTL_info, on = .(SOMAmer.ID, Sentinel.variant=`Sentinel.variant*`), type := `cis/.trans`]

# Filter to aptamers of interest
pQTLs[sinfo, on = .(SOMAmer.ID=SOMAMER_ID), Aptamer := Aptamer]
pQTLs <- pQTLs[Aptamer %in% prs_assocs$Aptamer]

# Filter columns
pQTLs <- pQTLs[, .(Aptamer, type, rsid=Conditional.variant, chr=Conditional.variant.Chr, pos=Conditional.variant.Pos)]

################################################################
# Get conditionally independent trans-pQTLs for aptamers not reported
# in Sun et al. 2018
#################################################################

not_mapped <- sinfo[!(trans_mapped), SOMAMER_ID]

new_trans <- foreach(this_somaid = not_mapped, .combine=rbind) %:% 
  foreach(this_chr = 1:22, .combine=rbind) %dopar% {
    this_seqid <- sinfo[SOMAMER_ID == this_somaid, Aptamer]
    this_cis_chr <- sinfo[Aptamer == this_seqid, chr]
    this_cis_tss <- as.integer(sinfo[Aptamer == this_seqid, start])
    this_target <- sinfo[Aptamer == this_seqid, Target]
    this_gene <- sinfo[Aptamer == this_seqid, Gene]
    this_uniprot <- sinfo[Aptamer == this_seqid, UniProt]

    pqtl_ss <- fread(sprintf("data/full_pQTL_summary_stats/%s/%s_chrom_%s_meta_1.tbl.gz", this_somaid, this_somaid, this_chr))
    pqtl_ss <- pqtl_ss[, .(Target = this_target, UniProt = this_uniprot, Gene = this_gene, Aptamer = this_seqid,
                           chr = chromosome, pos = position, EA = toupper(Allele1), OA = toupper(Allele2), Effect = Effect,
                           SE = StdErr, P = 10^(`log(P)`))]
    pqtl_ss <- pqtl_ss[!(chr == this_cis_chr & pos >= this_cis_tss - 1e6 & pos <= this_cis_tss + 1e6)]
    pqtl_ss[P < 1.5e-11]
}

fwrite(new_trans, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/new_trans_pQTLs.txt")

# Filter to top hit per locus
top_hits <- new_trans[,.SD[which.min(P)], by=.(Aptamer, chr)]

# Add top hits for aptamers of interest
pQTLs <- rbind(pQTLs, 
  top_hits[Aptamer %in% prs_assocs$Aptamer,.(Aptamer, type="trans", rsid=NA, chr, pos)]
)


################################################################
# Get top cis-pQTL for aptamers with no cis-pQTLs in Sun et al.
#################################################################

new_cis <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
new_cis[sinfo, on = .(SOMAMER_ID), Aptamer := Aptamer]
new_cis <- new_cis[Aptamer %in% prs_assocs$Aptamer]
new_cis <- new_cis[,.SD[which.min(P)],by=.(Aptamer)]
new_cis <- new_cis[,.(Aptamer, type="cis", rsid=NA, chr, pos)]
new_cis <- new_cis[!(Aptamer %in% pQTLs[type == "cis", Aptamer])]

pQTLs <- rbind(pQTLs, new_cis)

################################################################
# Extract dosages of pQTLs
#################################################################

dosage_info <- foreach(this_chr = pQTLs[,unique(chr)], .combine=rbind) %dopar% {
  # Get variant ids for extraction
  pQTL_pos <- pQTLs[chr == this_chr, pos]
  pfile <- sprintf("data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed/plink_format/pgen/impute_dedup_%s_interval", this_chr)
  pvar <- fread(sprintf("%s.pvar", pfile))
  pQTL_varid <- pvar[POS %in% pQTL_pos, .(ID)]

  # Write out variant list
  fwrite(pQTL_varid, col.names=FALSE, quote=FALSE, file=sprintf("analyses/pub/cardiometabolic_proteins/review2/pQTLs_chr%s.txt", this_chr))

  # Extract dosages
  cmd <- "plink2 --pfile"
  cmd[2] <- pfile
  cmd[3] <- sprintf("--extract analyses/pub/cardiometabolic_proteins/review2/pQTLs_chr%s.txt", this_chr)
  cmd[4] <- "--export A"
  cmd[5] <- sprintf("--out analyses/pub/cardiometabolic_proteins/review2/pQTLs_chr%s", this_chr)
  system(paste(cmd, collapse=" "), wait=TRUE)

  # Load dosage info
  di <- fread(sprintf("analyses/pub/cardiometabolic_proteins/review2/pQTLs_chr%s.raw", this_chr), nrows=0)
  di <- data.table(colid=names(di)[-(1:6)])
  di[, rsid := gsub("_.*", "", colid)]
  di[, dosage_allele := gsub(".*_", "", colid)]
  di[pvar[POS %in% pQTL_pos], on = .(rsid=ID), c("chr", "pos", "ALT", "REF") := .(this_chr, POS, ALT, REF)]
  return(di)
}

# Get allele frequencies
dosage_maf <- foreach(this_chr = dosage_info[,unique(chr)], .combine=rbind) %dopar% {
  this_pos <- dosage_info[chr == this_chr, pos]
  snpstats <- fread(sprintf("data/INTERVAL/reference_files/imputed_genotypes/impute_%s_interval.snpstats", this_chr))
  snpstats[position %in% this_pos, .(chr=this_chr, pos=position, minor_allele, major_allele, MAF)]
}
dosage_info[dosage_maf, on = .(chr, pos, dosage_allele = minor_allele), dosage_freq := MAF]
dosage_info[dosage_maf, on = .(chr, pos, dosage_allele = major_allele), dosage_freq := 1 - MAF]

# Add to pQTL info
pQTLs[, source := ifelse(is.na(rsid), "new", "Sun et al.")]
pQTLs[, rsid := NULL]
pQTLs <- pQTLs[dosage_info, on = .(chr, pos)]

# Load dosages
cmerge <- function(a1, a2) { merge(a1, a2, by="IID") }
dosages <- foreach(this_chr = pQTLs[,unique(chr)], .combine=cmerge) %dopar% {
  ds <- fread(sprintf("analyses/pub/cardiometabolic_proteins/review2/pQTLs_chr%s.raw", this_chr))
  ds[, c("FID", "PAT", "MAT", "SEX", "PHENOTYPE") := NULL]
  ds
}

# Clean up temporary files
foreach(this_chr = pQTLs[,unique(chr)]) %do% {
  system(sprintf("rm -f analyses/pub/cardiometabolic_proteins/review2/pQTLs_chr%s.txt", this_chr), wait=TRUE)
  system(sprintf("rm -f analyses/pub/cardiometabolic_proteins/review2/pQTLs_chr%s.log", this_chr), wait=TRUE)
  system(sprintf("rm -f analyses/pub/cardiometabolic_proteins/review2/pQTLs_chr%s.raw", this_chr), wait=TRUE)
}

# Write out combined dosage file
fwrite(dosages, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/prs_aptamer_QTL_dosages.txt")

# Filter to analysis subcohort
dosages <- dosages[!prev, on = .(IID)][IID %in% unique(soma$IID)]

# Compute dosage frequencies in this subgroup
dmaf2 <- data.table(colid=names(dosages)[-1], soma_freq =sapply(dosages[,-1], function(x) { sum(x) / (length(x)*2) }))
pQTLs[dmaf2, on = .(colid), soma_freq := soma_freq]

# Write out pQTL info file
fwrite(pQTLs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/pQTL_dosage_info.txt")

################################################################
# Test association between PRS and protein adjusting for pQTLs
#################################################################

# Load PRS
prs <- rbind(idcol="PRS",
  CAD_PRS = fread("analyses/GRS_profiles/CAD_metaGRS/profile.sscore.gz"),
  IS_PRS = fread("analyses/GRS_profiles/Stroke_metaGRS/profile.sscore.gz"),
  T2D_PRS = fread("analyses/GRS_profiles/T2D_2018/profile.sscore.gz"),
  Afib_PRS = fread("analyses/GRS_profiles/Afib_2018/profile.sscore.gz"),
  CKD_PRS = fread("analyses/GRS_profiles/CKD_2019/profile.sscore.gz")
)

# Load PCs
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

# Adjust PRS for PCs
prs <- prs[pcs, on = .(IID), nomatch=0]
prs[, prs_adj_pcs := lm(scale(score_sum) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +
                                           PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=PRS]
prs <- prs[, .(IID, PRS, prs_adj_pcs)]

# Orient dosages to their minor allele
pQTLs[, other_allele := ALT]
for (this_colid in pQTLs[soma_freq > 0.5, unique(colid)]) {
  dosages[[this_colid]] <- (dosages[[this_colid]] - 2)*-1
  pQTLs[colid == this_colid, c("dosage_allele", "other_allele", "dosage_freq", "soma_freq") := 
        .(other_allele, dosage_allele, 1 - dosage_freq, 1 - soma_freq)]
}

# Fit models adjusting for pQTLs
pqtl_adj <- foreach(test_idx = prs_assocs[,.I], .combine=rbind) %do% {
  # Extract info about this test
  this_prs <- prs_assocs[test_idx, PRS]
  this_target <- prs_assocs[test_idx, Target]
  this_uniprot <- prs_assocs[test_idx, UniProt]
  this_gene <- prs_assocs[test_idx, Gene]
  this_seqid <- prs_assocs[test_idx, Aptamer]
  this_aptvar <- sinfo[Aptamer == this_seqid, variable]

  if (!(this_seqid %in% pQTLs$Aptamer)) {
    return(NULL)
  }

  this_pqtls <- pQTLs[Aptamer == this_seqid, colid]
  this_rsid <- pQTLs[Aptamer == this_seqid, rsid]
  this_chr <- pQTLs[Aptamer == this_seqid, chr]
  this_pos <- pQTLs[Aptamer == this_seqid, pos]
  this_ea <- pQTLs[Aptamer == this_seqid, dosage_allele]
  this_oa <- pQTLs[Aptamer == this_seqid, other_allele]
  this_eaf <- pQTLs[Aptamer == this_seqid, soma_freq]
  this_type <- pQTLs[Aptamer == this_seqid, type]
  this_source <- pQTLs[Aptamer == this_seqid, source]

  # build model data.table
  dat <- prs[PRS == this_prs]
  dat <- dat[soma[variable == this_aptvar], on = .(IID), nomatch=0] 
  dat <- dat[dosages[,.SD,.SDcols=c("IID", this_pqtls)], on = .(IID), nomatch=0]

  # Build model formula:
  f <- sprintf("scale(soma_ivt_adj_batch) ~ scale(prs_adj_pcs) + %s", paste(this_pqtls, collapse=" + "))
  
  # Fit model
  l1 <- lm(as.formula(f), data=dat)
  cf <- coef(summary(l1))
  ci <- confint(l1)

  # Return coefficients for PRS and pQTLs
  data.table(PRS=this_prs, Target=this_target, UniProt=this_uniprot, Gene=this_gene, Aptamer=this_seqid,
             coefficient = c("PRS", this_rsid), chr = c("", this_chr), pos = c("", this_pos),
             effect_allele = c("", this_ea), other_allele = c("", this_oa),
             frequency = c("", this_eaf), type = c("", this_type), source = c("", this_source),
             Beta = cf[-1, 1], SE = cf[-1, 2], L95 = ci[-1, 1], U95 = ci[-1, 2], P = cf[-1, 4])
}

# Order rows for output
test_order <- pqtl_adj[coefficient == "PRS"][order(-abs(Beta))][order(P)][order(PRS)]
pqtl_adj <- pqtl_adj[order(-abs(Beta))][order(P)][order(source)][order(type)]
pqtl_adj <- pqtl_adj[test_order[,.(PRS, Aptamer)], on = .(PRS, Aptamer)]

fwrite(pqtl_adj, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/pqtl_adj_assocs.txt")

# Build comparison table for plotting
pqtl_adj_prot <- pqtl_adj[coefficient == "PRS", .(
  Beta = mean(Beta), SE = mean(SE), L95 = mean(L95), U95 = mean(U95), P = prot_pvalue(P, Beta)
), by = .(PRS, Target, UniProt, Gene)]
pqtl_adj_prot[, has_pQTLs := TRUE]

comp <- merge(prs_assocs, pqtl_adj_prot, by = c("PRS", "Target", "UniProt", "Gene"), all.x=TRUE,
              suffixes=c("", ".adj"))
comp[is.na(has_pQTLs), c("has_pQTLs", "Beta.adj", "L95.adj", "U95.adj") := .(FALSE, Beta, L95, U95)]
comp <- comp[,.(PRS, Target, UniProt, Gene, has_pQTLs, Beta, L95, U95, Beta.adj, L95.adj, U95.adj)]
comp <- unique(comp)

# Plot 
rdt <- comp[, .(x=c(min(c(L95, L95.adj)), 0,
                    max(c(U95, U95.adj)))*1.05,
                ymax=c(min(c(L95, L95.adj)), 0,
                       max(c(U95, U95.adj)))*1.05,
                ymin=c(0,0,0))]

g <- ggplot(comp) +
  aes(x = Beta, y = Beta.adj,
      xmin = L95, ymin = L95.adj,
      xmax = U95, ymax = U95.adj,
      color = has_pQTLs, fill = has_pQTLs) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
  geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.7, size=0.5) +
  geom_errorbar(width=0, alpha=0.7, size=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#bdbdbd") +
  geom_point(shape = 21, size=1.2, color="#00000000") +
  facet_wrap(~ PRS, nrow=1) +
  scale_x_continuous(name = "Beta (95% CI) for the PRS", expand=c(0,0)) +
  scale_y_continuous(name = "Beta (95% CI) adjusting for pQTLs") +
  scale_fill_manual(guide=FALSE,
    values=c("TRUE"="#1f78b4", "FALSE"="#525252")) +
  scale_color_manual(guide=FALSE,
    values=c("TRUE"="#a6cee3", "FALSE"="#969696")) +
  theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank()
  )
ggsave(g, width=7.2, height=2.4, file="analyses/pub/cardiometabolic_proteins/review2/pqtl_adj.pdf", useDingbats=FALSE)

