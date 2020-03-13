library(data.table)
library(foreach)

# Load somalogic IIDs for filtering
prot_iid <- fread("analyses/processed_traits/somalogic_proteins/traits.tsv")
prot_iid <- prot_iid[, unique(IID)]

# Load PCs so we can adjust when calculating correlations:
pcs <- fread("data/INTERVAL/reference_files/imputed_genotypes/annot_INT_50PCs_pcs.txt")
pcs <- pcs[, .(IID=ID, PC_1, PC_2, PC_3, PC_4, PC_5, PC_6, PC_7, PC_8, PC_9, PC_10)]

cor <- function(mat) {
  cmat <- matrix(NA, nrow=ncol(mat), ncol=ncol(mat), dimnames=list(colnames(mat), colnames(mat)))
  pmat <- cmat
  for (ci in colnames(mat)) {
    for (ri in colnames(mat)) {
      cmat[ri, ci] <- cor.test(mat[,ri], mat[,ci])$estimate
      pmat[ri, ci] <- cor.test(mat[,ri], mat[,ci])$p.value
    }
  }
  list(cmat, pmat)
}
        
# T2D PRS comparison
t2d <- rbind(
 "My_T2D"=fread("analyses/GRS_profiles/T2D_2018/profile.sscore.gz"),
 "Khera2018_T2D"=fread("analyses/GRS_profiles/Khera2018_T2D_PGS000014/profile.sscore.gz"),
 "Mahajan2018_T2D"=fread("analyses/GRS_profiles/Mahajan2018_T2D_PGS000036/profile.sscore.gz"),
 idcol="PRS")
t2d <- t2d[IID %in% prot_iid]
t2d[, score := scale(score_sum), by=PRS]
t2d <- t2d[pcs, on = .(IID), nomatch=0]
t2d[, adj_score := lm(score ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=PRS]
t2d <- dcast(t2d, IID ~ PRS, value.var="adj_score")
t2d <- as.matrix(t2d, rownames="IID")
cor(t2d)

# Afib PRS comparison
afib <- rbind(
 "My_afib"=fread("analyses/GRS_profiles/Afib_2018/profile.sscore.gz"),
 "Khera2018_Afib"=fread("analyses/GRS_profiles/Khera2018_afib_PGS000016/profile.sscore.gz"),
 "Weng2017_Afib"=fread("analyses/GRS_profiles/Weng2017_afib_PGS000035/profile.sscore.gz"),
 "Stroke_metaGRS_Afib"=fread("analyses/GRS_profiles/Afib_small/profile.sscore.gz"),
 idcol="PRS")
afib <- afib[IID %in% prot_iid]
afib[, score := scale(score_sum), by=PRS]
afib <- afib[pcs, on = .(IID), nomatch=0]
afib[, adj_score := lm(score ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=PRS]
afib <- dcast(afib, IID ~ PRS, value.var="adj_score")
afib <- as.matrix(afib, rownames="IID")
cor(afib)

# Stroke PRS comparison
stroke <- rbind(
 "My_stroke"=fread("analyses/GRS_profiles/StrokeAS_UKBv3_2018/profile.sscore.gz"),
 "Stroke_metaGRS"=fread("analyses/GRS_profiles/Stroke_metaGRS/profile.sscore.gz"),
 "Stroke_metaGRS_AnyStroke"=fread("analyses/GRS_profiles/StrokeAS_2018/profile.sscore.gz"),
 "RuttenJacobs2018_Stroke"=fread("analyses/GRS_profiles/RuttenJacobs2018_Stroke_PGSuncurated/profile.sscore.gz"),
 idcol="PRS")
stroke <- stroke[IID %in% prot_iid]
stroke[, score := scale(score_sum), by=PRS]
stroke <- stroke[pcs, on = .(IID), nomatch=0]
stroke[, adj_score := lm(score ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=PRS]
stroke <- dcast(stroke, IID ~ PRS, value.var="adj_score")
stroke <- as.matrix(stroke, rownames="IID")
cor(stroke)

# CKD PRS comparison
ckd <- rbind(
  "My_CKD"=fread("analyses/GRS_profiles/CKD_2019/profile.sscore.gz"),
  "Wuttke2019_eGFR"=fread("analyses/GRS_profiles/Wuttke2019_eGFR_PGSuncurated/profile.sscore.gz"),
  idcol="PRS")
ckd <- ckd[IID %in% prot_iid]
ckd[, score := scale(score_sum), by=PRS]
ckd <- ckd[pcs, on = .(IID), nomatch=0]
ckd[, adj_score := lm(score ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by=PRS]
ckd <- dcast(ckd, IID ~ PRS, value.var="adj_score")
ckd <- as.matrix(ckd, rownames="IID")
cor(ckd)


