library(data.table)
library(openxlsx)
source("src/utilities/prot_pval.R")
source("src/utilities/format_pval.R")

# Load and curate associations from our paper
hes_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")

# Extract protein level estimates
hes_assocs <- unique(hes_assocs[endpoint == "Diabetes",.(Protein=Gene, logHR=log(HR), logHR.L95=log(HR.L95),
                                logHR.U95=log(HR.U95), HR.P, HR.Nominal=HR.P<0.05, HR.Bonferroni=HR.P<0.05/42)])

# Filter to T2D PRS associated proteins
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[PRS == "T2D_PRS" & FDR < 0.05, .(Protein=Gene)])
hes_assocs <- hes_assocs[prs_assocs, on = .(Protein)]

# Load in prevalent and incident associations adjusted for age and sex from Gudmundsdottir et al. 2020
prev <- read.xlsx("data/Gudmundsdottir_et_al_2020/Gudmundsdottir_et_al_SupplementaryTables_030320.xlsx", sheet="Table S3", startRow=2)
setDT(prev)
prev <- prev[,c(1,3:6)]
setnames(prev, c("Protein", "AGES.Prev.Beta", "AGES.Prev.SE", "AGES.Prev.P", "AGES.Prev.Bonf"))
prev <- prev[AGES.Prev.Bonf < 0.05]
prev <- prev[,.(AGES.Prev.Beta=mean(AGES.Prev.Beta), AGES.Prev.SE=mean(AGES.Prev.SE), AGES.Prev.P=prot_pvalue(AGES.Prev.P, AGES.Prev.Beta)), by=Protein] 
prev[, AGES.Prev.L95 := AGES.Prev.Beta - qnorm(0.975) * AGES.Prev.SE]
prev[, AGES.Prev.U95 := AGES.Prev.Beta + qnorm(0.975) * AGES.Prev.SE]
prev <- prev[,.(Protein, AGES.Prev.Beta, AGES.Prev.SE, AGES.Prev.L95, AGES.Prev.U95, AGES.Prev.P)]

inci <- read.xlsx("data/Gudmundsdottir_et_al_2020/Gudmundsdottir_et_al_SupplementaryTables_030320.xlsx", sheet="Table S6", startRow=2)
setDT(inci)
inci <- inci[,c(1,3:6)]
setnames(inci, c("Protein", "AGES.Inci.Beta", "AGES.Inci.SE", "AGES.Inci.P", "AGES.Inci.Bonf"))
inci <- inci[AGES.Inci.Bonf < 0.05]
inci <- inci[,.(AGES.Inci.Beta=mean(AGES.Inci.Beta), AGES.Inci.SE=mean(AGES.Inci.SE), AGES.Inci.P=prot_pvalue(AGES.Inci.P, AGES.Inci.Beta)), by=Protein] 
inci[, AGES.Inci.L95 := AGES.Inci.Beta - qnorm(0.975) * AGES.Inci.SE]
inci[, AGES.Inci.U95 := AGES.Inci.Beta + qnorm(0.975) * AGES.Inci.SE]
inci <- inci[,.(Protein, AGES.Inci.Beta, AGES.Inci.SE, AGES.Inci.L95, AGES.Inci.U95, AGES.Inci.P)]

prev2 <- read.xlsx("data/Gudmundsdottir_et_al_2020/Gudmundsdottir_et_al_SupplementaryTables_030320.xlsx", sheet="Table S5", startRow=3)
setDT(prev2)
prev2 <- prev2[,c(1:5)]
setnames(prev2, c("Protein", "QMDiab.Prev.Beta", "QMDiab.Prev.SE", "QMDiab.Prev.P", "QMDiab.Prev.Bonf"))
prev2 <- prev2[QMDiab.Prev.Bonf < 0.05]
prev2 <- prev2[,.(QMDiab.Prev.Beta=mean(QMDiab.Prev.Beta), QMDiab.Prev.SE=mean(QMDiab.Prev.SE), QMDiab.Prev.P=prot_pvalue(QMDiab.Prev.P, QMDiab.Prev.Beta)), by=Protein] 
prev2[, QMDiab.Prev.L95 := QMDiab.Prev.Beta - qnorm(0.975) * QMDiab.Prev.SE]
prev2[, QMDiab.Prev.U95 := QMDiab.Prev.Beta + qnorm(0.975) * QMDiab.Prev.SE]
prev2 <- prev2[,.(Protein, QMDiab.Prev.Beta, QMDiab.Prev.SE, QMDiab.Prev.L95, QMDiab.Prev.U95, QMDiab.Prev.P)]

# Make comparison table and compare directional consistency
comp <- copy(hes_assocs)
comp <- merge(comp, prev, by = "Protein", all.x=TRUE)
comp[, AGES.Prev.Consistent := sign(AGES.Prev.Beta) == sign(logHR)]
comp <- merge(comp, prev2, by = "Protein", all.x=TRUE)
comp[, QMDiab.Prev.Consistent := sign(QMDiab.Prev.Beta) == sign(logHR)]
comp <- merge(comp, inci, by = "Protein", all.x=TRUE)
comp[, AGES.Inci.Consistent := sign(AGES.Inci.Beta) == sign(logHR)]

# Add mediation analysis
mediation <- fread("analyses/pub/cardiometabolic_proteins/review4/protein_mediation.txt")
mediation <- unique(mediation[PRS == "T2D_PRS" & Mediation == "natural indirect effect", 
  .(Protein=Gene, Mediation.PTE=PTE, Mediation.PTE.L95=PTE.L95, Mediation.PTE.U95=U95,
    Mediation.P=P, Mediation.Nominal=P<0.05, Mediation.Bonferroni=P<0.05/42)])
comp <- merge(comp, mediation, by="Protein")

# Add MR analysis
mr <- fread("analyses/pub/cardiometabolic_proteins/review2/mr_summary_all_proteins.txt")
comp[mr[PRS == "T2D_2018"], on = .(Protein=Gene), c("MR.logOR", "MR.L95", "MR.U95", "MR.P", "MR.Pleiotropy") :=
                                                  .(mr_estimate, mr_L95, mr_U95, mr_pval, pleiotropy_pval)]
comp[, MR.Nominal := MR.P < 0.05]
comp[, MR.Bonferroni := MR.P < 0.05/12]

# Add MR analysis from Gudmundsdottir et al.
gud_mr <- read.xlsx("data/Gudmundsdottir_et_al_2020/Gudmundsdottir_et_al_SupplementaryTables_030320.xlsx", sheet="Table S10", startRow=2)
setDT(gud_mr)
gud_mr <- gud_mr[,c(1,4:7, 26)]
setnames(gud_mr, c("Protein", "AGES.MR.Beta", "AGES.MR.SE", "AGES.MR.P", "AGES.MR.FDR", "AGES.MR.Consistent"))
gud_mr[, AGES.MR.L95 := AGES.MR.Beta - qnorm(0.975) * AGES.MR.SE]
gud_mr[, AGES.MR.U95 := AGES.MR.Beta + qnorm(0.975) * AGES.MR.SE]
gud_mr[, AGES.MR.Consistent := AGES.MR.Consistent == "Yes"]
gud_mr <- gud_mr[,.(Protein, AGES.MR.Beta, AGES.MR.SE, AGES.MR.L95, AGES.MR.U95, AGES.MR.P, AGES.MR.FDR, AGES.MR.Consistent)]

comp <- merge(comp, gud_mr, by="Protein", all.x=TRUE)

# Order rows
comp <- comp[order(HR.P)]

# write out
system("mkdir -p analyses/pub/cardiometabolic_proteins/review4", wait=TRUE)
fwrite(comp, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/gudmundsdottir_etal_T2D_comparison.txt")

# Make a plot comparing our HRs to the ORs in Gudmundsdottir et al.
gg_dt <- rbind(idcol="cohort",
  AGES.Prev=prev[,.(Protein, Beta=AGES.Prev.Beta, L95=AGES.Prev.L95, U95=AGES.Prev.U95)],
  QMDiab.Prev=prev2[,.(Protein, Beta=QMDiab.Prev.Beta, L95=QMDiab.Prev.L95, U95=QMDiab.Prev.U95)],
  AGES.Inci=inci[,.(Protein, Beta=AGES.Inci.Beta, L95=AGES.Inci.L95, U95=AGES.Inci.U95)]
)
gg_dt <- merge(hes_assocs[,.(Protein, logHR, L95=logHR.L95, U95=logHR.U95, HR.Nominal, HR.Bonferroni)], 
                gg_dt, by="Protein", suffixes=c("", ".Gudmundsdottir"))

gg_dt[, HR := exp(logHR)]
gg_dt[, L95 := exp(L95)]
gg_dt[, U95 := exp(U95)]
gg_dt[, OR.Gudmundsdottir := exp(Beta)]
gg_dt[, L95.Gudmundsdottir := exp(L95.Gudmundsdottir)]
gg_dt[, U95.Gudmundsdottir := exp(U95.Gudmundsdottir)]

gg_dt[, HR.type := fcase(HR.Bonferroni, "P < 0.001", !(HR.Bonferroni) & HR.Nominal, "P < 0.05", default = "P > 0.05")]
       

g <- ggplot(gg_dt) + 
   aes(y = HR, ymin = L95, ymax = U95,
       x = OR.Gudmundsdottir, xmin=L95.Gudmundsdottir, xmax=U95.Gudmundsdottir,
       shape = HR.Nominal, fill = HR.type) +
   geom_hline(yintercept=1, linetype=2) + geom_vline(xintercept=1, linetype=2) +
   geom_errorbarh(height=0) + geom_errorbar(width=0) + geom_point() +
   scale_shape_manual(values=c("TRUE"=23, "FALSE"=21)) +
   scale_fill_manual(values=c("P < 0.001"="yellow", "P < 0.05"="gray", "P > 0.05"="white")) +
   scale_x_log10("Odds Ratio (95% CI)", breaks=c(0.3, 0.5, 1, 2, 3)) +
   scale_y_log10("Hazard Ratio (95% CI)", breaks=c(0.3, 0.5, 1, 2, 3)) +
   facet_grid(~ cohort) +
   theme_bw() + 
   theme(legend.position="none", axis.text=element_text(size=6), axis.title=element_text(size=7),
         strip.text=element_text(size=7))
ggsave(g, width=7.2, height=2.8, file="analyses/pub/cardiometabolic_proteins/review4/gudmundsdottir_etal_T2D_comparison.pdf")

# Build contigency table for fisher exact testing:
overlaps <- rbind(use.names=FALSE,
  data.table(rowname="FDR < 0.05 with T2D PGS", row_total=comp[,.N], colname="Proteome-wide significant for prevalent T2D", col_total=prev[,.N], overlap=comp[AGES.Prev.P < 0.05, .N]),
  data.table("FDR < 0.05 with T2D PGS", comp[,.N], "Proteome-wide significant for incident T2D", inci[,.N], comp[AGES.Inci.P < 0.05, .N]),
  data.table("FDR < 0.05 with T2D PGS", comp[,.N], "Causal P < 0.05 in Mendelian randomization", gud_mr[,.N], comp[AGES.MR.P < 0.05, .N]),
  data.table("FDR < 0.05 with T2D PGS", comp[,.N], "Causal FDR < 0.05 in Mendelian randomization", gud_mr[AGES.MR.FDR < 0.05, .N], comp[AGES.MR.FDR < 0.05, .N]),

  data.table("P < 0.05 for incident T2D", comp[HR.P < 0.05,.N], "Proteome-wide significant for prevalent T2D", prev[,.N], comp[HR.P < 0.05 & AGES.Prev.P < 0.05, .N]),
  data.table("P < 0.05 for incident T2D", comp[HR.P < 0.05,.N], "Proteome-wide significant for incident T2D", inci[,.N], comp[HR.P < 0.05 & AGES.Inci.P < 0.05, .N]),
  data.table("P < 0.05 for incident T2D", comp[HR.P < 0.05,.N], "Causal P < 0.05 in Mendelian randomization", gud_mr[,.N], comp[HR.P < 0.05 & AGES.MR.P < 0.05, .N]),
  data.table("P < 0.05 for incident T2D", comp[HR.P < 0.05,.N], "Causal FDR < 0.05 in Mendelian randomization", gud_mr[AGES.MR.FDR < 0.05, .N], comp[HR.P < 0.05 & AGES.MR.FDR < 0.05, .N]),
   
  data.table("P < 0.05/42 for incident T2D", comp[HR.P < 0.05/42,.N], "Proteome-wide significant for prevalent T2D", prev[,.N], comp[HR.P < 0.05/42 & AGES.Prev.P < 0.05, .N]),
  data.table("P < 0.05/42 for incident T2D", comp[HR.P < 0.05/42,.N], "Proteome-wide significant for incident T2D", inci[,.N], comp[HR.P < 0.05/42 & AGES.Inci.P < 0.05, .N]),
  data.table("P < 0.05/42 for incident T2D", comp[HR.P < 0.05/42,.N], "Causal P < 0.05 in Mendelian randomization", gud_mr[,.N], comp[HR.P < 0.05/42 & AGES.MR.P < 0.05, .N]),
  data.table("P < 0.05/42 for incident T2D", comp[HR.P < 0.05/42,.N], "Causal FDR < 0.05 in Mendelian randomization", gud_mr[AGES.MR.FDR < 0.05, .N], comp[HR.P < 0.05/42 & AGES.MR.FDR < 0.05, .N]),
  
  data.table("Causal P < 0.05 in mediation analysis", comp[Mediation.P < 0.05,.N], "Proteome-wide significant for prevalent T2D", prev[,.N], comp[Mediation.P < 0.05 & AGES.Prev.P < 0.05, .N]),
  data.table("Causal P < 0.05 in mediation analysis", comp[Mediation.P < 0.05,.N], "Proteome-wide significant for incident T2D", inci[,.N], comp[Mediation.P < 0.05 & AGES.Inci.P < 0.05, .N]),
  data.table("Causal P < 0.05 in mediation analysis", comp[Mediation.P < 0.05,.N], "Causal P < 0.05 in Mendelian randomization", gud_mr[,.N], comp[Mediation.P < 0.05 & AGES.MR.P < 0.05, .N]),
  data.table("Causal P < 0.05 in mediation analysis", comp[Mediation.P < 0.05,.N], "Causal FDR < 0.05 in Mendelian randomization", gud_mr[AGES.MR.FDR < 0.05, .N], comp[Mediation.P < 0.05 & AGES.MR.FDR < 0.05, .N]),
   
  data.table("Causal P < 0.05/42 in mediation analysis", comp[Mediation.P < 0.05/42,.N], "Proteome-wide significant for prevalent T2D", prev[,.N], comp[Mediation.P < 0.05/42 & AGES.Prev.P < 0.05, .N]),
  data.table("Causal P < 0.05/42 in mediation analysis", comp[Mediation.P < 0.05/42,.N], "Proteome-wide significant for incident T2D", inci[,.N], comp[Mediation.P < 0.05/42 & AGES.Inci.P < 0.05, .N]),
  data.table("Causal P < 0.05/42 in mediation analysis", comp[Mediation.P < 0.05/42,.N], "Causal P < 0.05 in Mendelian randomization", gud_mr[,.N], comp[Mediation.P < 0.05/42 & AGES.MR.P < 0.05, .N]),
  data.table("Causal P < 0.05/42 in mediation analysis", comp[Mediation.P < 0.05/42,.N], "Causal FDR < 0.05 in Mendelian randomization", gud_mr[AGES.MR.FDR < 0.05, .N], comp[Mediation.P < 0.05/42 & AGES.MR.FDR < 0.05, .N]),

  data.table("Causal P < 0.05 in Mendelian randomization analysis", comp[MR.P < 0.05,.N], "Proteome-wide significant for prevalent T2D", prev[,.N], comp[MR.P < 0.05 & AGES.Prev.P < 0.05, .N]),
  data.table("Causal P < 0.05 in Mendelian randomization analysis", comp[MR.P < 0.05,.N], "Proteome-wide significant for incident T2D", inci[,.N], comp[MR.P < 0.05 & AGES.Inci.P < 0.05, .N]),
  data.table("Causal P < 0.05 in Mendelian randomization analysis", comp[MR.P < 0.05,.N], "Causal P < 0.05 in Mendelian randomization", gud_mr[,.N], comp[MR.P < 0.05 & AGES.MR.P < 0.05, .N]),
  data.table("Causal P < 0.05 in Mendelian randomization analysis", comp[MR.P < 0.05,.N], "Causal FDR < 0.05 in Mendelian randomization", gud_mr[AGES.MR.FDR < 0.05, .N], comp[MR.P < 0.05 & AGES.MR.FDR < 0.05, .N]),
   
  data.table("Causal P < 0.05/12 in Mendelian randomization analysis", comp[MR.P < 0.05/12,.N], "Proteome-wide significant for prevalent T2D", prev[,.N], comp[MR.P < 0.05/12 & AGES.Prev.P < 0.05, .N]),
  data.table("Causal P < 0.05/12 in Mendelian randomization analysis", comp[MR.P < 0.05/12,.N], "Proteome-wide significant for incident T2D", inci[,.N], comp[MR.P < 0.05/12 & AGES.Inci.P < 0.05, .N]),
  data.table("Causal P < 0.05/12 in Mendelian randomization analysis", comp[MR.P < 0.05/12,.N], "Causal P < 0.05 in Mendelian randomization", gud_mr[,.N], comp[MR.P < 0.05/12 & AGES.MR.P < 0.05, .N]),
  data.table("Causal P < 0.05/12 in Mendelian randomization analysis", comp[MR.P < 0.05/12,.N], "Causal FDR < 0.05 in Mendelian randomization", gud_mr[AGES.MR.FDR < 0.05, .N], comp[MR.P < 0.05/12 & AGES.MR.FDR < 0.05, .N])
)

# Gudmundsdottir use SomaLogic V4 platform, using Entrez gene IDs for distinct proteins.
total_overlap <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
total_overlap <- total_overlap[Type == "Protein" & update_anno == "V3_in_V4", length(unique(Gene.Name))]
overlaps[, total := total_overlap]

# get cells for fisher test
overlaps[, true_true := overlap]
overlaps[, true_false := col_total - overlap]
overlaps[, false_true := row_total - overlap]
overlaps[, false_false := total - true_true - true_false - false_true]

# Fisher exact test
overlaps[, Fisher.Test.P := fisher.test(matrix(c(true_true, false_true, true_false, false_false), nrow=2), alternative="greater")$p.value, by=.(rowname, colname)]

# Output
overlaps[, text := sprintf("%s (P = %s)", overlap, my_format_pval(Fisher.Test.P))]

contingency <- dcast(overlaps, rowname ~ colname, value.var="text")
contingency <- as.matrix(contingency, rownames=1)
contingency <- rbind(contingency, unique(overlaps[, .(colname, col_total)])[colnames(contingency), on = .(colname)][["col_total"]])
contingency <- cbind(contingency, c(unique(overlaps[, .(rowname, row_total)])[rownames(contingency), on = .(rowname), nomatch=0][["row_total"]], total_overlap))

rownames(contingency)[nrow(contingency)] <- "Column totals"
colnames(contingency)[ncol(contingency)] <- "Row totals"

contingency <- contingency[
  c("FDR < 0.05 with T2D PGS", "P < 0.05 for incident T2D", "P < 0.05/42 for incident T2D", "Causal P < 0.05 in mediation analysis", 
    "Causal P < 0.05/42 in mediation analysis", "Causal P < 0.05 in Mendelian randomization analysis", "Causal P < 0.05/12 in Mendelian randomization analysis", "Column totals"),
  c("Proteome-wide significant for prevalent T2D", "Proteome-wide significant for incident T2D", "Causal P < 0.05 in Mendelian randomization",
    "Causal FDR < 0.05 in Mendelian randomization", "Row totals")
]
contingency <- as.data.table(contingency, keep.rownames="rownames")

fwrite(contingency, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review4/gudmundsdottir_etal_T2D_comparison_fisher_test.txt")

