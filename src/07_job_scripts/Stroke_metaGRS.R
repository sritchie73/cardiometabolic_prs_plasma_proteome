library(data.table)

out_dir <- "analyses/mendelian_randomisation/Stroke_metaGRS/gwas_summary_stats"
ss_dir <- "data/GWAS_summary_statistics/MEGASTROKE"
name <- "StrokeIS"

dir.create(sprintf("%s/%s/", out_dir, name), recursive=TRUE, showWarnings=FALSE)

ss <- fread(sprintf("%s/MEGASTROKE.2.IS.EUR.GC_filtered_X_nocases_het.TBL.gz", ss_dir))
info <- data.table(samples=ss[,max(TotalSampleSize)], cases=ss[,max(TotalCases)], 
                   controls=ss[, max(TotalSampleSize) - max(TotalCases)])
fwrite(info, sep="\t", quote=FALSE, file=sprintf("%s/%s/info.txt", out_dir, name))

ss <- ss[, .(var_id=MarkerName, EA=toupper(Allele1), OA=toupper(Allele2),
              MAF=Freq1, effect=Effect, se=StdErr, P=`P-value`)]
setkey(ss, var_id, EA, OA)

for (chr in 1:22) {
  bim <- fread(sprintf("analyses/processed_genotypes/impute_chr%s_interval_filtered.bim", chr))
  setnames(bim, c("chr", "var_id", "cm", "pos", "a1", "a2"))
  setkey(bim, var_id, a1, a2)

  chr_ss <- rbind(
		ss[bim, on = .(var_id, EA=a1, OA=a2), nomatch=0],
		ss[bim, on = .(var_id, EA=a2, OA=a1), nomatch=0])
  chr_ss <- chr_ss[, .(var_id, chr, pos, EA, OA, MAF, effect, se, P)]
 
  fwrite(chr_ss, file=sprintf("%s/%s/chr%s.txt", out_dir, name, chr), sep="\t", quote=FALSE)
  system(command=sprintf("gzip %s/%s/chr%s.txt", out_dir, name, chr), wait=TRUE)
}

