library(data.table)

# Create new environment for temporary processing
metabo.tmp <- new.env()
with(metabo.tmp, {

# Load raw metabolon data
raw <- rbind(idcol="batch", fill=TRUE, use.names=TRUE,
   fread("data/INTERVAL/post_qc_data/phenotype/metabolon_metabolomics/qc/B1_ScaledImpData_with_missingness.txt"),
   fread("data/INTERVAL/post_qc_data/phenotype/metabolon_metabolomics/qc/B2_ScaledImpData_with_missingness.txt"))

# Load in GWAS QC data and filter to samples and metabolites passing QC:
gwasqc <- fread("data/INTERVAL/post_qc_data/phenotype/metabolon_metabolomics/gwasqc/metabolon_gwasqc.csv")
setnames(raw, "INTERVAL_ALIQUOT_ID", "aliquot_id")
metab_cols <- names(raw)[grepl("^[0-9]*$", names(raw))]
setnames(raw, metab_cols, sprintf("m%05d", as.integer(metab_cols)))
raw <- raw[,.SD,.SDcols=c(names(gwasqc), "PLATE_ID", "batch")]
raw <- raw[aliquot_id %in% gwasqc$aliquot_id]

# Filter to creatinine measure
raw <- raw[,.(aliquot_id, creatinine=m00513, PLATE_ID, batch)]
raw <- raw[!is.na(creatinine)]

# Load phenotype data and obtain covariates
pheno <- fread("data/INTERVAL/project_1074/INTERVALdata_17JUL2018.csv")
idmap <- fread("data/INTERVAL/project_1074/omicsMap.csv")
pheno[idmap, on = .(identifier), Metabolon_met_gwasQC_bl := Metabolon_met_gwasQC_bl]
pheno[, daysToProcess_bl := as.IDate(processDate_bl, format="%d%b%Y") - as.IDate(attendanceDate, format="%d%b%Y")]
raw[pheno, on = .(aliquot_id = Metabolon_met_gwasQC_bl), c("daysToProcess", "age", "sex") := .(daysToProcess_bl, agePulse, sexPulse)]

# Drop withdrawn samples:
raw <- raw[aliquot_id %in% pheno$Metabolon_met_gwasQC_bl]

# Map to genetic identifiers
raw[idmap, on = .(aliquot_id=Metabolon_met_gwasQC_bl), IID := Affymetrix_gwasQC_bl]
raw <- raw[!is.na(IID)]

# Log transform creatinine for adjustment
raw[, log_creatinine := log(creatinine)]

# Adjust creatinine for technical covariates:
raw[, log_adjusted := lm(log_creatinine ~ factor(PLATE_ID) + factor(batch) + daysToProcess)$residuals]

# Put back onto raw units - residuals always have mean (and sum) = 0, and SD is always 
# smaller than SD of the original data (the difference in SD before and after adjustment,
# known as the residual standard error in econometrics, is proportional to the variance 
# explained by the dependent variables; in this case the technical factors we're adjusting for.
# Provided we trust the original measurements are well calibrated (in terms of population mean)
# we can rescale the residuals to absolute units by shifting the mean back to its original value,
# then reverse the log transformation
raw[, log_adjusted := log_adjusted + mean(log_creatinine)]
raw[, adjusted := exp(log_adjusted)]

# Filter columns
raw <- raw[,.(IID, creatinine=adjusted, age, sex)]

# Compute eGFR using the CKD-EPI equation
raw[sex == 1 & creatinine <= 0.9, eGFR := 141 * (creatinine/0.9)^-0.411 * 0.993^floor(age)]
raw[sex == 1 & creatinine > 0.9, eGFR := 141 * (creatinine/0.9)^-1.209 * 0.993^floor(age)]
raw[sex == 2 & creatinine <= 0.7, eGFR := 141 * (creatinine/0.9)^-0.329 * 0.993^floor(age)]
raw[sex == 2 & creatinine > 0.7, eGFR := 141 * (creatinine/0.9)^-1.209 * 0.993^floor(age)]

# Compute CKD stage - not likely relevant as it may misclassify elderly people (who have reduced kidney function)
# as having disease, and there is a strong age trend to the below
raw[eGFR >= 90, CKD := "Normal"] # CKD stage 1 if evidence of kidney damage
raw[eGFR < 90 & eGFR >= 60, CKD := "Reduced function"] # Requires evidence of Kidney Damage to also be present for CKD stage 2
raw[eGFR < 60 & eGFR >= 30, CKD := "Moderate"] 
raw[eGFR < 30 & eGFR >= 15, CKD := "Severe"]
raw[eGFR < 15, CKD := "Kidney failure"]

# Finalise
})
crea <- metabo.tmp$raw[,.(IID, creatinine, eGFR, CKD)]
rm(metabo.tmp)

