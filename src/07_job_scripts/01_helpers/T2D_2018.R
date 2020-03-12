library(data.table)

out_dir <- "analyses/mendelian_randomisation/T2D_2018/gwas_summary_stats"
ss_dir <- "data/GWAS_summary_statistics/Type2Diabetes_2018"
name <- "T2D"

dir.create(sprintf("%s/%s", out_dir, name), recursive=TRUE, showWarnings=FALSE)

ss <- fread(cmd=sprintf("unzip -p %s/Mahajan.NatGenet2018b.T2D.European.zip", ss_dir))
ss <- ss[, .(chr=Chr, pos=Pos, EA, OA=NEA, MAF=EAF, effect=Beta, se=SE, P=Pvalue)]
setkey(ss, chr, pos, EA, OA)

for (chr in 1:22) {
  bim <- fread(sprintf("analyses/processed_genotypes/impute_chr%s_interval_filtered.bim", chr))
  setnames(bim, c("chr", "var_id", "cm", "pos", "a1", "a2"))
  setkey(bim, chr, pos, a1, a2)

  chr_ss <- rbind(
    ss[bim, on = .(chr, pos, EA=a1, OA=a2), nomatch=0],
    ss[bim, on = .(chr, pos, EA=a2, OA=a1), nomatch=0])
  chr_ss <- chr_ss[, .(var_id, chr, pos, EA, OA, MAF, effect, se, P)]

  fwrite(chr_ss, file=sprintf("%s/%s/chr%s.txt", out_dir, name, chr), sep="\t", quote=FALSE)
  system(command=sprintf("gzip %s/%s/chr%s.txt", out_dir, name, chr), wait=TRUE)
}

info <- data.table(samples=898130, cases=74124, controls=824006)
fwrite(info, sep="\t", quote=FALSE, file=sprintf("%s/%s/info.txt", out_dir, name))

name <- "T2DadjBMI"
dir.create(sprintf("%s/%s", out_dir, name), recursive=TRUE, showWarnings=FALSE)

ss <- fread(cmd=sprintf("unzip -p %s/Mahajan.NatGenet2018b.T2Dbmiadj.European.zip", ss_dir))
ss <- ss[, .(chr=Chr, pos=Pos, EA, OA=NEA, MAF=EAF, effect=Beta, se=SE, P=Pvalue)]
setkey(ss, chr, pos, EA, OA)

for (chr in 1:22) {
  bim <- fread(sprintf("analyses/processed_genotypes/impute_chr%s_interval_filtered.bim", chr))
  setnames(bim, c("chr", "var_id", "cm", "pos", "a1", "a2"))
  setkey(bim, chr, pos, a1, a2)

  chr_ss <- rbind(
    ss[bim, on = .(chr, pos, EA=a1, OA=a2), nomatch=0],
    ss[bim, on = .(chr, pos, EA=a2, OA=a1), nomatch=0])
  chr_ss <- chr_ss[, .(var_id, chr, pos, EA, OA, MAF, effect, se, P)]

  fwrite(chr_ss, file=sprintf("%s/%s/chr%s.txt", out_dir, name, chr), sep="\t", quote=FALSE)
  system(command=sprintf("gzip %s/%s/chr%s.txt", out_dir, name, chr), wait=TRUE)
}

info <- data.table(samples=898130, cases=74124, controls=824006)
fwrite(info, sep="\t", quote=FALSE, file=sprintf("%s/%s/info.txt", out_dir, name))
