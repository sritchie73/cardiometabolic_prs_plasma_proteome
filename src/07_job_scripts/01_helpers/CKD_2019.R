library(data.table)

out_dir <- "analyses/mendelian_randomisation/CKD_2019/gwas_summary_stats"
ss_dir <- "data/GWAS_summary_statistics/CKDGen_2019"
name <- "CKD"

dir.create(sprintf("%s/%s/", out_dir, name), recursive=TRUE, showWarnings=FALSE)

ss <- fread(sprintf("%s/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz", ss_dir))
ss <- ss[, .(var_id=RSID, chr=Chr, pos=Pos_b37,
              EA=toupper(Allele1), OA=toupper(Allele2),
              MAF=Freq1, effect=Effect, se=StdErr, P=`P-value`)]
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

info <- data.table(samples=480698, cases=41395, controls=439303)
fwrite(info, sep="\t", quote=FALSE, file=sprintf("%s/%s/info.txt", out_dir, name))



