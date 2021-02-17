library(data.table)
library(openxlsx)
library(foreach)
library(doMC)
library(ggplot2)
library(ggrastr)
library(pheatmap)
library(RColorBrewer)

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

# Load PRS
prs <- rbind(idcol="PRS",
  CAD_PRS = fread("analyses/GRS_profiles/CAD_metaGRS/profile.sscore.gz"),
  IS_PRS = fread("analyses/GRS_profiles/Stroke_metaGRS/profile.sscore.gz"),
  T2D_PRS = fread("analyses/GRS_profiles/T2D_2018/profile.sscore.gz"),
  AF_PRS = fread("analyses/GRS_profiles/Afib_2018/profile.sscore.gz"),
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

# Load SomaLogic aptamer levels
soma <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")

# Adjust for batch
batch <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")
soma <- soma[batch, on = .(IID), nomatch=0]
soma[, value := lm(value ~ factor(batch))$residuals, by=variable]
soma[, soma_ivt_adj_batch := scale(value)]
soma <- soma[,.(IID, variable, soma_ivt_adj_batch)]

# Filter prs table to samples with somalogic measures
prs <- prs[IID %in% unique(soma$IID)]
prs[, prs_adj_pcs := scale(prs_adj_pcs)] # standardise to SD = 1, mean = 0

# Drop prevalent cases
prs <- prs[!prev, on = .(IID)]
soma <- soma[!prev, on = .(IID)]

# And withdrawn samples
idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
idmap <- fread(idmap_file)

prs <- prs[IID %in% idmap$Affymetrix_gwasQC_bl]
soma <- soma[IID %in% idmap$Affymetrix_gwasQC_bl]

# Load protein information.
sinfo <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

# Filter to aptamers passing QC
sinfo <- sinfo[Type == "Protein"]

# Select columns
sinfo <- sinfo[Type == "Protein", .(variable, SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, 
                                    Gene=Gene.Name, Entrez_id=Entrez.Gene.ID, chr, start, end, Cross_Reactivity=Characterization.Info,
                                    Mass_Spec_Confirmation=Mass.Spec.Confirmation.in.Matrix)]

# Curate aptamer sensitivity and specificity information
emilsson <- rbind(fill=TRUE,
  as.data.table(read.xlsx("data/Emilsson_etal_2018/supp_tables.xlsx", sheet="Table S3", startRow=2)),
  as.data.table(read.xlsx("data/Emilsson_etal_2018/supp_tables.xlsx", sheet="Table S4", startRow=2)))
soma_by_gene <- sinfo[,.(Gene=strsplit(Gene, "\\|")[[1]]),by=variable]
soma_by_gene[emilsson, on = .(Gene=Gene.Symbol), Mass_Spec_Confirmation := TRUE]
soma_by_gene <- soma_by_gene[!is.na(Mass_Spec_Confirmation)]

sinfo[Mass_Spec_Confirmation != "", Mass_Spec_Confirmation := "SomaLogic"]
sinfo[soma_by_gene, on = .(variable), Mass_Spec_Confirmation := ifelse(Mass_Spec_Confirmation == "", "Emilsson et al. 2018", Mass_Spec_Confirmation)]

# Are the aptamers supported by cis pQTLs?
cis_pQTLs <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
cis_pQTLs <- unique(cis_pQTLs[,.(SOMAMER_ID)])
sinfo[, cis_pQTL := ""]
sinfo[cis_pQTLs, on = .(SOMAMER_ID), cis_pQTL := "yes"]

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

# Filter to aptamers passing QC
soma <- soma[variable %in% sinfo$variable]

####################################
# Get PRS to protein associations
####################################

# Fit associations for each aptamer and PRS.
apt_assocs <- foreach(this_prs = prs[,unique(PRS)], .combine=rbind) %:% 
  foreach(this_apt = soma[,unique(variable)], .combine=rbind) %dopar% {
    dat <- merge(soma[variable == this_apt], prs[PRS == this_prs], by="IID")
    l1 <- lm(soma_ivt_adj_batch ~ prs_adj_pcs, data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)
    data.table(PRS = this_prs, Aptvar = this_apt, Beta=cf[2,1], SE=cf[2,2], L95=ci[2,1], U95=ci[2,2], P=cf[2,4])
}

# Add protein information
apt_assocs <- apt_assocs[sinfo, on = .(Aptvar=variable)]

# Aggregate at the protein target level
prot_assocs <- apt_assocs[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)),
                          by = .(PRS, Target, UniProt, Gene, Entrez_id, chr, start, end)]
prot_assocs[, FDR := p.adjust(P, method="fdr"), by=.(PRS)]

# Combine with aptamer associations
assocs <- merge(prot_assocs, apt_assocs, by=c("PRS", "Target", "UniProt", "Gene", "Entrez_id", "chr", "start", "end"),
                suffixes=c("", ".Aptamer"))

# Row and column order
assocs <- assocs[order(-abs(Beta.Aptamer))][order(P.Aptamer)][order(-abs(Beta))][order(P)][order(PRS)]
assocs <- assocs[,.(PRS, Target, UniProt, Gene, Entrez_id, chr, start, end, Beta, SE, L95, U95, P, FDR, 
                    Aptamer, Beta.Aptamer, SE.Aptamer, L95.Aptamer, U95.Aptamer, P.Aptamer, Cross_Reactivity,
                    Mass_Spec_Confirmation, cis_pQTL)]

# Write out supplementary table
fwrite(assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")

########################################
# Create QQ plots
########################################
qq_assocs <- prot_assocs[,.(PRS, Target, UniProt, Gene, P, FDR)]
qq_assocs[, sig := "FALSE"]
qq_assocs[FDR < 0.1, sig := "Suggestive"]
qq_assocs[FDR < 0.05, sig := "Significant"]
qq_assocs[, PRS := gsub("_", " ", PRS)]

# Generate expected distribution of p-values
qq_assocs <- qq_assocs[order(P)][order(PRS)]
qq_assocs[, expected := ppoints(.N), by = PRS]

# Transform for plot
qq_assocs[, expected := -log10(expected)]
qq_assocs[, observed := -log10(P)]

# Get top 5
anno <- qq_assocs[, .SD[1:5], by=PRS]

# Plot
g <- ggplot(qq_assocs) +
  aes(x=expected, y=observed) +
  geom_abline(intercept=0, slope=1, linetype=2, color="#737373") +
  geom_line(color="black", size=0.5) +
  geom_point_rast(shape=19, aes(color=sig), size=0.8, raster.width=1.270, raster.height=1.281, raster.dpi=1200) +
  geom_hline(color="#0571b0", linetype=2, yintercept=-log10(0.05)) +
  facet_wrap(~ PRS, nrow=1) +
  scale_x_continuous(name = "Expected -log10 p-values") +
  scale_y_continuous(name = "Observed -log10 p-values") +
  scale_color_manual(values=c("FALSE"="#000000", "Suggestive"="#4575b4",
                              "Significant"="#d73027"), guide=FALSE) +
  theme_bw() + theme(
    strip.text=element_text(size=8), axis.title=element_text(size=10),
    axis.text=element_text(size=8), panel.grid=element_blank()
  )
ggsave(g, width=7.2, height=2, file="analyses/pub/cardiometabolic_proteins/review2/prs_to_prot_assoc_qq.pdf", useDingbats=FALSE)

############################################
# Heatmap of PRS to protein associations
############################################

# Build data.table for heatmap. Show all proteins where there
# was FDR < 0.05 for any PRS.
heatmap_prots <- prot_assocs[FDR < 0.05, unique(Gene)]
heatmap_dt <- prot_assocs[Gene %in% heatmap_prots]
heatmap_dt[,PRS := gsub("_", " ", PRS)]

# Create wide tables for heatmap
beta_mat <- as.matrix(dcast(heatmap_dt, PRS ~ Gene, value.var="Beta"), rownames="PRS")
pval_mat <- as.matrix(dcast(heatmap_dt, PRS ~ Gene, value.var="P"), rownames="PRS")
fdr_mat <- as.matrix(dcast(heatmap_dt, PRS ~ Gene, value.var="FDR"), rownames="PRS")

# Create text matrix
text_mat <- matrix("", nrow=nrow(beta_mat), ncol=ncol(beta_mat),
                      dimnames=dimnames(beta_mat))
text_mat[pval_mat < 0.05] <- "o"
text_mat[fdr_mat < 0.1] <- "*"
text_mat[fdr_mat < 0.05] <- "#"

# Order heatmap
prs_order <- rownames(beta_mat)[hclust(dist(beta_mat))$order]
prot_order <- foreach(pIdx = seq_along(prs_order), .combine=c) %do% {
  # First show positive the negative associations, keeping aside any proteins
  # associated with the next PRS. Drop any proteins associated with the previous
  # prs in the order so we don't have duplicates.
  this_prs <- heatmap_dt[PRS == prs_order[pIdx] & FDR < 0.05]
  if (pIdx > 1) {
    assoc_with_last <- heatmap_dt[PRS == prs_order[pIdx - 1] & FDR < 0.05, .(Gene)]
    this_prs <- this_prs[!assoc_with_last, on = .(Gene)]
  }
  if (pIdx != length(prs_order)) {
    assoc_with_next <- heatmap_dt[PRS == prs_order[pIdx + 1] & FDR < 0.05, .(Gene)]
    both_prs <- this_prs[assoc_with_next, on = .(Gene), nomatch = 0]
    this_prs <- this_prs[!assoc_with_next, on = .(Gene)]
  }
  po <- c(this_prs[Beta > 0][order(-Beta), Gene], this_prs[Beta < 0][order(Beta), Gene])
  if (pIdx != length(prs_order)) {
    po <- c(po, both_prs[Beta > 0][order(-Beta), Gene], both_prs[Beta < 0][order(Beta), Gene])
  }
  return(po)
}
beta_mat <- beta_mat[prs_order, prot_order]
text_mat <- text_mat[prs_order, prot_order]

# Split heatmap into positive and negative associations
pos_prot <- heatmap_dt[FDR < 0.05 & Beta > 0, Gene]
pos_prot <- intersect(prot_order, pos_prot)
neg_prot <- heatmap_dt[FDR < 0.05 & Beta < 0, Gene]
neg_prot <- intersect(prot_order, neg_prot)

pal <- rev(brewer.pal(name="RdYlBu", n=8))
pal <- c(pal[1:4], "#FFFFFF", pal[5:8])

pheatmap(beta_mat[, pos_prot], display_numbers=text_mat[, pos_prot],
  cluster_cols=FALSE, cluster_rows=FALSE,
  color=colorRampPalette(pal)(255), breaks=seq(-0.12, 0.12, length=256),
  legend_breaks=seq(-0.12, 0.10, by=0.02), border_color=NA,
  fontsize_row=8, fontsize_col=8, fontsize_number=8,
  cellwidth=9.5, cellheight=10,
  filename="analyses/pub/cardiometabolic_proteins/review2/sig_heatmap_pos.pdf")

pheatmap(beta_mat[, neg_prot], display_numbers=text_mat[, neg_prot],
  cluster_cols=FALSE, cluster_rows=FALSE,
  color=colorRampPalette(pal)(255), breaks=seq(-0.12, 0.12, length=256),
  legend_breaks=seq(-0.12, 0.10, by=0.02), border_color=NA,
  fontsize_row=8, fontsize_col=8, fontsize_number=8,
  cellwidth=9.5, cellheight=10,
  filename="analyses/pub/cardiometabolic_proteins/review2/sig_heatmap_neg.pdf")

############################################
# Test multivariable associations with SHBG
############################################

# Test CAD, IS, and T2D PRS for joint association with SHBG
dat <- dcast(prs, IID ~ PRS, value.var="prs_adj_pcs")
dat <- dat[soma[variable %like% "shbg"], on = .(IID)]
shbg_assocs <- foreach(apt_var = unique(dat$variable), .combine=rbind) %do% {
  l1 <- lm(soma_ivt_adj_batch ~ CAD_PRS + IS_PRS + T2D_PRS, data=dat[variable == apt_var])
  cf <- coef(summary(l1))
  ci <- confint(l1)
  data.table(Aptamer=apt_var, PRS=c("CAD_PRS", "IS_PRS", "T2D_PRS"), Beta=cf[-1,1], 
             SE=cf[-1,2], L95=ci[-1,1], U95=ci[-1,2], P=cf[-1,4])
}
shbg_assocs[sinfo, on = .(Aptamer=variable), Aptamer := i.Aptamer]
shbg_assocs_agg <- shbg_assocs[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)), by=PRS]
shbg_assocs <- merge(shbg_assocs_agg, shbg_assocs, by="PRS", suffixes=c("", ".Aptamer"))

fwrite(shbg_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/shbg_multi_prs_assocs.txt")

############################################
# Sensitivity analysis
############################################

# Load phenotype data and filter to analysis cohort
pheno <- fread("analyses/processed_traits/phenotypes.tsv")
pheno <- pheno[IID %in% unique(soma$IID)]

# Compute BMI
pheno[wt_bl == 777, wt_bl := NA] # bad coding
pheno[, bmi_raw := wt_bl/ht_bl^2]
pheno[, bmi := bmi_raw]
pheno[ht_bl < 1.47, bmi := NA_real_] # clinical cutoff for dwarfism
pheno[ht_bl > 2.1, bmi := NA_real_] # clinical cutoff for gigantism
pheno[wt_bl < 50 | wt_bl > 160, bmi := NA_real_] # NHS restrictions for weight

# Select columns
pheno <- pheno[,.(IID, bmi, attendanceDate, appointmentTime)]

# Convert dates and times to correct formats
pheno[, attendanceDate := as.IDate(attendanceDate, format="%d%B%Y")]
pheno[, appointmentTime := as.ITime(appointmentTime)]

# Bin into 10 bins of equal duration, with the largest bin in terms
# of sample size used as the reference group. If NAs are present, these
# are treated as an 11th bin (and cannot be the reference group).

bin_duration <- function(x, n=10) {
  # Determine duration cut points for bins
  min_x <- min(x, na.rm=TRUE)
  max_x <- max(x, na.rm=TRUE)
  bin_cuts <- seq(min_x, max_x, length=11)

  reclass <- if ("ITime" %in% class(x)) {
    as.ITime 
  } else if ("IDate" %in% class(x)) {
    as.IDate
  }

  # Compute number of samples in each bin
  dt <- data.table(x=x)
  bin_cuts <- data.table(bin_start=reclass(bin_cuts[1:n]), bin_end=reclass(bin_cuts[(1:n)+1]))
  bin_cuts[dt, on = .(bin_start <= x, bin_end > x), samples := .N, by=.(bin_start, bin_end)]
  bin_cuts[n, samples := samples + dt[x == max_x, .N]]

  # Assign bin number based on number of samples
  bin_cuts <- bin_cuts[order(-samples)]
  bin_cuts[, bin_number := .I]

  # If there are missing values, set as n+1th bin
  if (any(is.na(x))) {
    bin_cuts <- rbind(bin_cuts, data.table(samples=sum(is.na(x)), bin_number=n+1), fill=TRUE)
  }
 
  # Return bins
  bin_cuts
}

time_bins <- bin_duration(pheno$appointmentTime)
date_bins <- bin_duration(pheno$attendanceDate)

fwrite(time_bins, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/somalogic_samples_time_bins.txt")
fwrite(date_bins, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/somalogic_samples_date_bins.txt")

pheno[time_bins, on = .(appointmentTime >= bin_start, appointmentTime < bin_end), time_bin := bin_number]
pheno[time_bins[which.max(bin_end)], on = .(appointmentTime = bin_end), time_bin := bin_number]
pheno[is.na(appointmentTime), time_bin := 11]

pheno[date_bins, on = .(attendanceDate >= bin_start, attendanceDate < bin_end), date_bin := bin_number]
pheno[date_bins[which.max(bin_end)], on = .(attendanceDate = bin_end), date_bin := bin_number]

######################################################################
# Fit associations adjusting protein levels for time of blood draw
#####################################################################

# Adjust protein levels for time of day
soma[pheno, on = .(IID), time_bin := i.time_bin]
soma[, time_adj := scale(lm(soma_ivt_adj_batch ~ factor(time_bin))$residuals), by=variable]

# Fit associations for each aptamer and PRS.
time_apt_assocs <- foreach(this_prs = prs[,unique(PRS)], .combine=rbind) %:% 
  foreach(this_apt = soma[,unique(variable)], .combine=rbind) %dopar% {
    dat <- merge(soma[variable == this_apt], prs[PRS == this_prs], by="IID")
    l1 <- lm(time_adj ~ prs_adj_pcs, data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)
    data.table(PRS = this_prs, Aptvar = this_apt, Beta=cf[2,1], SE=cf[2,2], L95=ci[2,1], U95=ci[2,2], P=cf[2,4])
}

# Add protein information
time_apt_assocs <- time_apt_assocs[sinfo, on = .(Aptvar=variable)]

# Aggregate at the protein target level
time_prot_assocs <- time_apt_assocs[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)),
                                     by = .(PRS, Target, UniProt, Gene, Entrez_id, chr, start, end)]

# Combine with aptamer associations
time_assocs <- merge(time_prot_assocs, time_apt_assocs, by=c("PRS", "Target", "UniProt", "Gene", "Entrez_id", "chr", "start", "end"),
                suffixes=c("", ".Aptamer"))

# Row and column order
time_assocs <- time_assocs[order(-abs(Beta.Aptamer))][order(P.Aptamer)][order(-abs(Beta))][order(P)][order(PRS)]
time_assocs <- time_assocs[,.(PRS, Target, UniProt, Gene, Entrez_id, chr, start, end, Beta, SE, L95, U95, P, 
                    Aptamer, Beta.Aptamer, SE.Aptamer, L95.Aptamer, U95.Aptamer, P.Aptamer)]

# Write out supplementary table
fwrite(time_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_time_assocs.txt")

######################################################################
# Fit associations adjusting protein levels for date of blood draw
#####################################################################

# Adjust protein levels for date of day
soma[pheno, on = .(IID), date_bin := i.date_bin]
soma[, date_adj := scale(lm(soma_ivt_adj_batch ~ factor(date_bin))$residuals), by=variable]

# Fit associations for each aptamer and PRS.
date_apt_assocs <- foreach(this_prs = prs[,unique(PRS)], .combine=rbind) %:% 
  foreach(this_apt = soma[,unique(variable)], .combine=rbind) %dopar% {
    dat <- merge(soma[variable == this_apt], prs[PRS == this_prs], by="IID")
    l1 <- lm(date_adj ~ prs_adj_pcs, data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)
    data.table(PRS = this_prs, Aptvar = this_apt, Beta=cf[2,1], SE=cf[2,2], L95=ci[2,1], U95=ci[2,2], P=cf[2,4])
}

# Add protein information
date_apt_assocs <- date_apt_assocs[sinfo, on = .(Aptvar=variable)]

# Aggregate at the protein target level
date_prot_assocs <- date_apt_assocs[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)),
                                     by = .(PRS, Target, UniProt, Gene, Entrez_id, chr, start, end)]

# Combine with aptamer associations
date_assocs <- merge(date_prot_assocs, date_apt_assocs, by=c("PRS", "Target", "UniProt", "Gene", "Entrez_id", "chr", "start", "end"),
                suffixes=c("", ".Aptamer"))

# Row and column order
date_assocs <- date_assocs[order(-abs(Beta.Aptamer))][order(P.Aptamer)][order(-abs(Beta))][order(P)][order(PRS)]
date_assocs <- date_assocs[,.(PRS, Target, UniProt, Gene, Entrez_id, chr, start, end, Beta, SE, L95, U95, P, 
                    Aptamer, Beta.Aptamer, SE.Aptamer, L95.Aptamer, U95.Aptamer, P.Aptamer)]

# Write out supplementary table
fwrite(date_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_date_assocs.txt")

########################################################################################################
# Fit associations adjusting for BMI (may affect protein levels, or PRS association may be through BMI)
########################################################################################################

# Fit associations for each aptamer ~ PRS + BMI
bmi_apt_assocs <- foreach(this_prs = prs[,unique(PRS)], .combine=rbind) %:% 
  foreach(this_apt = soma[,unique(variable)], .combine=rbind) %dopar% {
    dat <- merge(soma[variable == this_apt], prs[PRS == this_prs], by="IID")
    dat[pheno, on = .(IID), bmi := i.bmi]
    dat <- dat[!is.na(bmi)]
    l1 <- lm(soma_ivt_adj_batch ~ prs_adj_pcs + scale(log(bmi)), data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)
    data.table(PRS = this_prs, Aptvar = this_apt, Coef=c("PRS", "BMI"), 
               Beta=cf[-1,1], SE=cf[-1,2], L95=ci[-1,1], U95=ci[-1,2], P=cf[-1,4])
}

# Add protein information
bmi_apt_assocs <- bmi_apt_assocs[sinfo, on = .(Aptvar=variable)]

# Aggregate at the protein target level
bmi_prot_assocs <- bmi_apt_assocs[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)),
                                     by = .(PRS, Coef, Target, UniProt, Gene, Entrez_id, chr, start, end)]

# Combine with aptamer associations
bmi_assocs <- merge(bmi_prot_assocs, bmi_apt_assocs, by=c("PRS", "Coef", "Target", "UniProt", "Gene", "Entrez_id", "chr", "start", "end"),
                suffixes=c("", ".Aptamer"))

# Row and column order
bmi_assocs <- bmi_assocs[order(-abs(Beta.Aptamer))][order(P.Aptamer)][order(-abs(Beta))][order(P)][order(Coef)][order(PRS)]
bmi_assocs <- bmi_assocs[,.(PRS, Coef, Target, UniProt, Gene, Entrez_id, chr, start, end, Beta, SE, L95, U95, P, 
                    Aptamer, Beta.Aptamer, SE.Aptamer, L95.Aptamer, U95.Aptamer, P.Aptamer)]

# Write out supplementary table
fwrite(bmi_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_bmi_assocs.txt")

########################################################################################################
# Fit associations between BMI on PRS and protein 
########################################################################################################

# Fit reverse associations for BMI ~ aptamer + PRS
bmi_apt_assocs2 <- foreach(this_prs = prs[,unique(PRS)], .combine=rbind) %:% 
  foreach(this_apt = soma[,unique(variable)], .combine=rbind) %dopar% {
    dat <- merge(soma[variable == this_apt], prs[PRS == this_prs], by="IID")
    dat[pheno, on = .(IID), bmi := i.bmi]
    dat <- dat[!is.na(bmi)]
    l1 <- lm(scale(log(bmi)) ~ prs_adj_pcs + soma_ivt_adj_batch, data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)
    data.table(PRS = this_prs, Aptvar = this_apt, Coef=c("PRS", "Aptamer"), 
               Beta=cf[-1,1], SE=cf[-1,2], L95=ci[-1,1], U95=ci[-1,2], P=cf[-1,4])
}

# Add protein information
bmi_apt_assocs2 <- bmi_apt_assocs2[sinfo, on = .(Aptvar=variable)]

# Aggregate at the protein target level
bmi_prot_assocs2 <- bmi_apt_assocs2[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)),
                                     by = .(PRS, Coef, Target, UniProt, Gene, Entrez_id, chr, start, end)]

# Combine with aptamer associations
bmi_assocs2 <- merge(bmi_prot_assocs2, bmi_apt_assocs2, by=c("PRS", "Coef", "Target", "UniProt", "Gene", "Entrez_id", "chr", "start", "end"),
                suffixes=c("", ".Aptamer"))

# Row and column order
bmi_assocs2 <- bmi_assocs2[order(-abs(Beta.Aptamer))][order(P.Aptamer)][order(-abs(Beta))][order(P)][order(Coef)][order(PRS)]
bmi_assocs2 <- bmi_assocs2[,.(PRS, Coef, Target, UniProt, Gene, Entrez_id, chr, start, end, Beta, SE, L95, U95, P, 
                    Aptamer, Beta.Aptamer, SE.Aptamer, L95.Aptamer, U95.Aptamer, P.Aptamer)]

# Write out supplementary table
fwrite(bmi_assocs2, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/bmi_on_prs_and_prot_assocs.txt")

######################################################
# Create comparison table
######################################################

# Get significant PRS to protein associations
bmi_comp <- prot_assocs[FDR < 0.05]
bmi_comp <- bmi_comp[,.(PRS, Target, UniProt, Gene, Beta, SE, L95, U95, P, FDR)]

# Add in BMI adjusted results
bmi_comp[bmi_prot_assocs[Coef == "PRS"], on = .(PRS, Target, UniProt, Gene), 
  paste0(c("Beta", "SE", "L95", "U95", "P"), ".BMIadj.PRS") := 
  .(i.Beta, i.SE, i.L95, i.U95, i.P)]

# Filter to associations attenuated by BMI adjustment
bmi_comp <- bmi_comp[P.BMIadj.PRS > 0.05]

# Add in BMI coefficient estimates
bmi_comp[bmi_prot_assocs[Coef == "BMI"], on = .(PRS, Target, UniProt, Gene), 
  paste0(c("Beta", "SE", "L95", "U95", "P"), ".BMIadj.BMI") := 
  .(i.Beta, i.SE, i.L95, i.U95, i.P)]

# Add in PRS and Protein coefficients for when BMI is the response
bmi_comp[bmi_prot_assocs2[Coef == "PRS"], on = .(PRS, Target, UniProt, Gene), 
  paste0("BMI.", c("Beta", "SE", "L95", "U95", "P"), ".PRS") := 
  .(i.Beta, i.SE, i.L95, i.U95, i.P)]
bmi_comp[bmi_prot_assocs2[Coef == "Aptamer"], on = .(PRS, Target, UniProt, Gene), 
  paste0("BMI.", c("Beta", "SE", "L95", "U95", "P"), ".Protein") := 
  .(i.Beta, i.SE, i.L95, i.U95, i.P)]

fwrite(bmi_comp, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/bmi_attenuated_assocs.txt")

######################################################
# Associations including prevalent cases
######################################################

# Reload somalogic data
soma <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")

# Adjust for batch
batch <- fread("analyses/processed_traits/somalogic_proteins/covariates.tsv")
soma <- soma[batch, on = .(IID), nomatch=0]
soma[, value := lm(value ~ factor(batch))$residuals, by=variable]
soma[, soma_ivt_adj_batch := scale(value)]
soma <- soma[,.(IID, variable, soma_ivt_adj_batch)]

# Drop withdrawn sample
soma <- soma[IID %in% idmap$Affymetrix_gwasQC_bl]

# Filter to aptamers passing QC
soma <- soma[variable %in% sinfo$variable]

# Fit associations for each aptamer and PRS.
prev_apt_assocs <- foreach(this_prs = prs[,unique(PRS)], .combine=rbind) %:% 
  foreach(this_apt = soma[,unique(variable)], .combine=rbind) %dopar% {
    dat <- merge(soma[variable == this_apt], prs[PRS == this_prs], by="IID")
    l1 <- lm(soma_ivt_adj_batch ~ prs_adj_pcs, data=dat)
    cf <- coef(summary(l1))
    ci <- confint(l1)
    data.table(PRS = this_prs, Aptvar = this_apt, Beta=cf[2,1], SE=cf[2,2], L95=ci[2,1], U95=ci[2,2], P=cf[2,4])
}

# Add protein information
prev_apt_assocs <- prev_apt_assocs[sinfo, on = .(Aptvar=variable)]

# Aggregate at the protein target level
prev_prot_assocs <- prev_apt_assocs[, .(Beta=mean(Beta), SE=mean(SE), L95=mean(L95), U95=mean(U95), P=prot_pvalue(P, Beta)),
                                     by = .(PRS, Target, UniProt, Gene, Entrez_id, chr, start, end)]
prev_prot_assocs[, FDR := p.adjust(P, method="fdr"), by=.(PRS)]

# Combine with aptamer associations
prev_assocs <- merge(prev_prot_assocs, prev_apt_assocs, by=c("PRS", "Target", "UniProt", "Gene", "Entrez_id", "chr", "start", "end"),
                     suffixes=c("", ".Aptamer"))

# Row and column order
prev_assocs <- prev_assocs[order(-abs(Beta.Aptamer))][order(P.Aptamer)][order(-abs(Beta))][order(P)][order(PRS)]
prev_assocs <- prev_assocs[,.(PRS, Target, UniProt, Gene, Entrez_id, chr, start, end, Beta, SE, L95, U95, P, FDR, 
															Aptamer, Beta.Aptamer, SE.Aptamer, L95.Aptamer, U95.Aptamer, P.Aptamer, Cross_Reactivity,
															Mass_Spec_Confirmation, cis_pQTL)]

# Write out supplementary table
fwrite(prev_assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs_with_prev.txt")

######################################################
# Create sensitivity analysis plot
######################################################

# Extract significant PRS to protein associations
sig_assocs <- prot_assocs[FDR < 0.05]

# Extract and combine sensitivity analyses
sens_dt <- rbind(idcol="covariate", fill=TRUE,
  "time_of_day"=time_prot_assocs[sig_assocs[, .(PRS, Target, UniProt, Gene)], on = .(PRS, Target, UniProt, Gene)],
  "date_of_donation"=date_prot_assocs[sig_assocs[, .(PRS, Target, UniProt, Gene)], on = .(PRS, Target, UniProt, Gene)],
  "bmi"=bmi_prot_assocs[sig_assocs[, .(PRS, Target, UniProt, Gene)], on = .(PRS, Target, UniProt, Gene)][Coef == "PRS"],
  "prevalent_disease"=prev_prot_assocs[sig_assocs[, .(PRS, Target, UniProt, Gene)], on = .(PRS, Target, UniProt, Gene)]
)

# Filter to minimal set of needed columns
sens_dt <- sens_dt[, .(covariate, PRS, Target, Beta, L95, U95)]
sig_assocs <- sig_assocs[,.(PRS, Target, Beta, L95, U95)]

# Build comparison table
sens_comp <- merge(sig_assocs, sens_dt, by=c("PRS", "Target"), suffixes=c("", ".adj"))

# Data table for ribbon
rdt <- sens_comp[, .(x=c(min(c(L95, L95.adj)), 0,
                         max(c(U95, U95.adj)))*1.05,
												 ymax=c(min(c(L95, L95.adj)), 0,
												  			max(c(U95, U95.adj)))*1.05,
												 ymin=c(0,0,0))]
# Plot
g <- ggplot(sens_comp) + 
  aes(x=Beta, xmin=L95, xmax=U95, y=Beta.adj, ymin=L95.adj, ymax=U95.adj, 
      color=covariate, fill=covariate) +
  geom_ribbon(data=rdt, aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="#f0f0f0") +
  geom_hline(yintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_vline(xintercept=0, color="#bdbdbd", linetype="dotted") +
  geom_errorbarh(height=0, alpha=0.7, size=0.5) +
  geom_errorbar(width=0, alpha=0.7, size=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#737373") +
  geom_point(shape = 21, size=1, color="#00000000") +
  facet_grid(covariate ~ PRS) +
  scale_x_continuous(name = "Beta (95% CI)", expand=c(0,0)) +
  scale_y_continuous(name = "Beta (95% CI)") +
  scale_color_manual(guide=FALSE, values=c("bmi"="#fdbb84", "date_of_donation"="#addd8e", "time_of_day"="#41b6c4", "prevalent_disease"="#807dba")) +
  scale_fill_manual(guide=FALSE, values=c("bmi"="#b30000", "date_of_donation"="#238b45", "time_of_day"="#225ea8", "prevalent_disease"="#54278f")) +
  theme_bw() + theme(
    axis.title=element_text(size=10), axis.text=element_text(size=8),
    panel.grid=element_blank()
  )
ggsave(g, width=7.2, height=7.2, useDingbats=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/sensitivity.pdf")
  
  
  
























  



