library(data.table)

# Load all sets of associations
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
polygenicity <- fread("analyses/pub/cardiometabolic_proteins/review2/polygenicity.txt")
hes_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")
mediation <- fread("analyses/pub/cardiometabolic_proteins/review3/protein_mediation.txt")

# Extract protein level estimates
prs_assocs <- unique(prs_assocs[,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])
hes_assocs <- unique(hes_assocs[,.(Endpoint=endpoint, Target, UniProt, Gene, HR, HR.L95, HR.U95, HR.P, HR.FDR)])
mediation <- unique(mediation[, .(PRS, Target, UniProt, Gene, Endpoint, Mediation, logOR, L95, U95, P, PTE)])

# Map between tables
map <- data.table(
  PRS = c("AF_PRS", "CAD_PRS", "CKD_PRS", "IS_PRS", "T2D_PRS"),
  Endpoint = c("Atrial fibrillation", "Myocardial infarction", "End stage renal disease", "Ischaemic stroke", "Diabetes")
)
prs_assocs[map, on = .(PRS), Endpoint := i.Endpoint]
hes_assocs[map, on = .(Endpoint), PRS := i.PRS]

# Format numbers
format_num <- function(x, centre=0) {
  x <- x - centre
  ifelse(abs(x) < 0.1, 
      sprintf("%.03f", round(x * 1000)/1000 + centre),
      sprintf("%.02f", round(x * 100)/100 + centre))
}

format_ci <- function(l95, u95, centre=0) {
  sprintf("[%s, %s]", format_num(l95, centre), format_num(u95, centre))
}

format_p <- function(p) {
  sapply(p, function(y) {
		fcase(
			y >= 0.1, sprintf("%.02f", round(y*100)/100),
			y >= 0.001, sprintf("%.03f", round(y*1000)/1000),
			default = gsub("e-0?", "x10-", format(y, digits=1, scientific=TRUE)))
  })
}

format_pct <- function(x) {
  ifelse(x < 0.1, 
      sprintf("%.01f%%", round(x * 1000)/10),
      sprintf("%d%%", round(x * 100)))

}

# Build table
comp <- prs_assocs[PRS.FDR < 0.05, .(PRS, Endpoint, Gene)]
comp[prs_assocs, on = .(PRS, Gene), PRS.Beta := format_num(i.PRS.Beta)]
comp[prs_assocs, on = .(PRS, Gene), PRS.95CI := format_ci(PRS.L95, PRS.U95)]
comp[prs_assocs, on = .(PRS, Gene), PRS.P := format_p(i.PRS.P)]
comp[prs_assocs, on = .(PRS, Gene), PRS.FDR := format_p(i.PRS.FDR)]
comp[polygenicity, on = .(PRS, Gene), Polygenicity := format_pct(pct_removed)]
comp[hes_assocs, on = .(Endpoint, Gene), HR := format_num(i.HR, 1)]
comp[hes_assocs, on = .(Endpoint, Gene), HR.95CI := format_ci(HR.L95, HR.U95, 1)]
comp[hes_assocs, on = .(Endpoint, Gene), HR.P := format_p(i.HR.P)]

med <- mediation[Mediation == "natural indirect effect"]
comp[med, on = .(PRS, Endpoint, Gene), Causal.OR := format_num(exp(logOR), 1)]
comp[med, on = .(PRS, Endpoint, Gene), Causal.95CI := format_ci(exp(L95), exp(U95), 1)]
comp[med, on = .(PRS, Endpoint, Gene), Causal.P := format_p(P)]
comp[med, on = .(PRS, Endpoint, Gene), Causal.PTE := format_pct(PTE)]

# Get row order:
c1 <- med[P < 0.05][order(P), .(PRS, Gene, Endpoint)]
c2 <- hes_assocs[comp[, .(PRS, Gene, Endpoint)], on = .(PRS, Gene, Endpoint), nomatch=0]
c2 <- c2[!c1, on = .(PRS, Gene, Endpoint)]
c2 <- c2[order(HR.P), .(PRS, Gene, Endpoint)]
c3 <- prs_assocs[PRS.FDR < 0.05][order(PRS.FDR), .(PRS, Gene, Endpoint)]
c3 <- c3[!c2, on = .(PRS, Gene, Endpoint)][!c1, on = .(PRS, Gene, Endpoint)]
ro <- rbind(c1, c2, c3)[order(PRS)]
comp <- comp[ro, on = .(PRS, Gene, Endpoint)]

fwrite(comp, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/compare_estimates.txt")

# Now do all the aptamer specific estimates
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
hes_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_biomarker_assocs.txt")
mediation <- fread("analyses/pub/cardiometabolic_proteins/review3/protein_mediation.txt")
bmi_mediation <- fread("analyses/pub/cardiometabolic_proteins/review3/bmi_mediation.txt")

# Extract aptamer level estimates
prs_assocs <- prs_assocs[,.(PRS, Gene, Aptamer, PRS.Beta=Beta.Aptamer, PRS.L95=L95.Aptamer, PRS.U95=U95.Aptamer, PRS.P=P.Aptamer, PRS.FDR=FDR)]
hes_assocs <- hes_assocs[,.(Endpoint=endpoint, Gene, Aptamer, HR=HR.aptamer, HR.L95=HR.L95.aptamer, HR.U95=HR.U95.aptamer, HR.P=HR.P.aptamer)]
mediation <- mediation[, .(PRS, Gene, Aptamer, Endpoint, Mediation, logOR=logOR.aptamer, L95=L95.aptamer, U95=U95.aptamer, P=P.aptamer, PTE=PTE.aptamer)]
bmi_mediation <- bmi_mediation[, .(PRS, Gene, Aptamer, Mediation, Beta=Beta.aptamer, L95=L95.aptamer, U95=U95.aptamer, P=P.aptamer, PTE=PTE.aptamer)]

# Build table
comp <- prs_assocs[PRS.FDR < 0.05, .(PRS, Gene, Aptamer)]
comp[map, on = .(PRS), Endpoint := Endpoint]
comp <- comp[,.(Endpoint, PRS, Gene, Aptamer)]
multapt <- comp[,.N,by=.(Endpoint, PRS, Gene)][N > 1, .(Endpoint, PRS, Gene)]
comp <- comp[multapt, on = .(Endpoint, PRS, Gene)]

comp[prs_assocs, on = .(PRS, Gene, Aptamer), PRS.Beta := format_num(i.PRS.Beta)]
comp[prs_assocs, on = .(PRS, Gene, Aptamer), PRS.95CI := format_ci(PRS.L95, PRS.U95)]
comp[prs_assocs, on = .(PRS, Gene, Aptamer), PRS.P := format_p(i.PRS.P)]

comp[bmi_mediation[Mediation == "natural indirect effect"], on = .(PRS, Gene, Aptamer), BMI.med.Beta := format_num(i.Beta)]
comp[bmi_mediation[Mediation == "natural indirect effect"], on = .(PRS, Gene, Aptamer), BMI.med.95CI := format_ci(L95, U95)]
comp[bmi_mediation[Mediation == "natural indirect effect"], on = .(PRS, Gene, Aptamer), BMI.med.P := format_p(P)]
comp[bmi_mediation[Mediation == "natural indirect effect"], on = .(PRS, Gene, Aptamer), BMI.med.PTE := format_pct(PTE)]

comp[bmi_mediation[Mediation == "natural direct effect"], on = .(PRS, Gene, Aptamer), BMI.ind.Beta := format_num(i.Beta)]
comp[bmi_mediation[Mediation == "natural direct effect"], on = .(PRS, Gene, Aptamer), BMI.ind.95CI := format_ci(L95, U95)]
comp[bmi_mediation[Mediation == "natural direct effect"], on = .(PRS, Gene, Aptamer), BMI.ind.P := format_p(P)]

comp[hes_assocs, on = .(Endpoint, Gene, Aptamer), HR := format_num(i.HR, 1)]
comp[hes_assocs, on = .(Endpoint, Gene, Aptamer), HR.95CI := format_ci(HR.L95, HR.U95, 1)]
comp[hes_assocs, on = .(Endpoint, Gene, Aptamer), HR.P := format_p(i.HR.P)]

med <- mediation[Mediation == "natural indirect effect"]
comp[med, on = .(PRS, Endpoint, Gene, Aptamer), Causal.OR := format_num(exp(logOR), 1)]
comp[med, on = .(PRS, Endpoint, Gene, Aptamer), Causal.95CI := format_ci(exp(L95), exp(U95), 1)]
comp[med, on = .(PRS, Endpoint, Gene, Aptamer), Causal.P := format_p(P)]
comp[med, on = .(PRS, Endpoint, Gene, Aptamer), Causal.PTE := format_pct(PTE)]

# Get row order:
comp <- comp[order(Aptamer)][order(Gene)][order(PRS)]
fwrite(comp, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/compare_aptamer_estimates.txt")

