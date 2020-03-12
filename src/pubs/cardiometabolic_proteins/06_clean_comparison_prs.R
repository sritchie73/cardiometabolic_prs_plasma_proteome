library(data.table)
library(foreach)
library(openxlsx)
source("src/utilities/flip_strand.R")

# Load bim files in case some columns are missing from curated PRSs
bim <- foreach(chrIdx = 1:22, .combine=rbind) %do% {
  fread(sprintf("analyses/processed_genotypes/impute_chr%s_interval_filtered.bim", chrIdx))
}
setnames(bim, c("chr", "rsid", "cm", "pos", "A1", "A2"))
setkey(bim, rsid, chr, pos, A1, A2)

# Convert each PRS to the common format for calculating levels in INTERVAL
dt <- fread("data/PGS_catalog/Khera2018_afib_PGS000016.txt.gz")
setnames(dt, c("chr", "pos", "effect_allele", "other_allele", "weight"))
dt[, rsid := sprintf("%s_%s_%s_%s", chr, pos, effect_allele, other_allele)]
dt <- dt[, .(rsid, chr, pos, effect_allele, other_allele, weight)]
dir.create("data/GRS_resources/Khera2018_afib_PGS000016")
fwrite(dt, sep="\t", quote=FALSE, file="data/GRS_resources/Khera2018_afib_PGS000016/grs_weights.txt")
system("gzip data/GRS_resources/Khera2018_afib_PGS000016/grs_weights.txt", wait=TRUE)

dt <- fread("data/PGS_catalog/Khera2018_T2D_PGS000014.txt.gz")
setnames(dt, c("chr", "pos", "effect_allele", "other_allele", "weight"))
dt[, rsid := sprintf("%s_%s_%s_%s", chr, pos, effect_allele, other_allele)]
dt <- dt[, .(rsid, chr, pos, effect_allele, other_allele, weight)]
dir.create("data/GRS_resources/Khera2018_T2D_PGS000014")
fwrite(dt, sep="\t", quote=FALSE, file="data/GRS_resources/Khera2018_T2D_PGS000014/grs_weights.txt")
system("gzip data/GRS_resources/Khera2018_T2D_PGS000014/grs_weights.txt", wait=TRUE)

dt <- fread("data/PGS_catalog/Mahajan2018_T2D_PGS000036.txt.gz")
setnames(dt, c("rsid", "chr", "pos", "effect_allele", "other_allele", "weight", "OR"))
dt <- dt[, .(rsid, chr, pos, effect_allele, other_allele, weight)]
dir.create("data/GRS_resources/Mahajan2018_T2D_PGS000036")
fwrite(dt, sep="\t", quote=FALSE, file="data/GRS_resources/Mahajan2018_T2D_PGS000036/grs_weights.txt")
system("gzip data/GRS_resources/Mahajan2018_T2D_PGS000036/grs_weights.txt", wait=TRUE)

dt <- fread("data/PGS_catalog/Weng2017_afib_PGS000035.txt.gz")
setnames(dt, c("rsid", "effect_allele", "weight"))
setkey(dt, rsid, effect_allele)
dt[bim, on = .(rsid, effect_allele=A1), c("chr", "pos", "other_allele") := .(i.chr, i.pos, A2)]
dt[bim, on = .(rsid, effect_allele=A2), c("chr", "pos", "other_allele") := .(i.chr, i.pos, A1)]
dt <- dt[order(pos)][order(chr)]
dt <- dt[, .(rsid, chr, pos, effect_allele, other_allele, weight)]
miss <- dt[is.na(chr), .N]
dt <- dt[!is.na(chr)]
dir.create("data/GRS_resources/Weng2017_afib_PGS000035")
fwrite(dt, sep="\t", quote=FALSE, file="data/GRS_resources/Weng2017_afib_PGS000035/grs_weights.txt")
system("gzip data/GRS_resources/Weng2017_afib_PGS000035/grs_weights.txt", wait=TRUE)
cat(miss, "\n", file="data/GRS_resources/Weng2017_afib_PGS000035/variants_not_passing_qc.txt")

dt <- fread("data/PGS_catalog/RuttenJacobs2018_Stroke_PGSuncarated.txt.gz")
setnames(dt, c("rsid", "effect_allele", "weight", "pval"))
setkey(dt, rsid, effect_allele)
dt[bim, on = .(rsid, effect_allele=A1), c("chr", "pos", "other_allele") := .(i.chr, i.pos, A2)]
dt[bim, on = .(rsid, effect_allele=A2), c("chr", "pos", "other_allele") := .(i.chr, i.pos, A1)]
dt <- dt[order(pos)][order(chr)]
dt <- dt[, .(rsid, chr, pos, effect_allele, other_allele, weight)]
miss <- dt[is.na(chr), .N]
dt <- dt[!is.na(chr)]
dir.create("data/GRS_resources/RuttenJacobs2018_Stroke_PGSuncurated")
fwrite(dt, sep="\t", quote=FALSE, file="data/GRS_resources/RuttenJacobs2018_Stroke_PGSuncurated/grs_weights.txt")
system("gzip data/GRS_resources/RuttenJacobs2018_Stroke_PGSuncurated/grs_weights.txt", wait=TRUE)
cat(miss, "\n", file="data/GRS_resources/RuttenJacobs2018_Stroke_PGSuncurated/variants_not_passing_qc.txt")

dt <- read.xlsx("data/PGS_catalog/Wuttke2019_supp_tables.xlsx", sheet="ST3", startRow=2)
dt <- as.data.table(dt)
dt <- dt[-309] # tailing comment
dt <- dt[`1-sided.p-value.MVP` < 0.05] # replicating eGFR SNPs
dt[,chr := gsub(":.*", "", `Chr/Pos.(b37)`)]
dt[,pos := gsub(".*:", "", `Chr/Pos.(b37)`)]
dt <- dt[,.(rsid=RS.number, chr, pos, effect_allele=EA, other_allele=NEA, weight=Effect)]
dt[, rsid := gsub(" .*$", "", rsid)]
dt[, effect_allele := toupper(effect_allele)]
dt[, other_allele := toupper(other_allele)]
dir.create("data/GRS_resources/Wuttke2019_eGFR_PGSuncurated")
fwrite(dt, sep="\t", quote=FALSE, file="data/GRS_resources/Wuttke2019_eGFR_PGSuncurated/grs_weights.txt")
system("gzip data/GRS_resources/Wuttke2019_eGFR_PGSuncurated/grs_weights.txt", wait=TRUE)









