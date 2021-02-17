library(data.table)

out_dir <- "analyses/mendelian_randomisation/Afib_2018/gwas_summary_stats"
ss_dir <- "data/GWAS_summary_statistics/Atrial_fib/"
name <- "Afib"

dir.create(sprintf("%s/%s/", out_dir, name), recursive=TRUE, showWarnings=FALSE)

ss <- fread(sprintf("%s/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz", ss_dir))
ss <- ss[,.(var_id=rs_dbSNP147, chr=as.numeric(CHR), pos=POS_GRCh37,
            EA=A2, OA=A1, MAF=Freq_A2, effect=Effect_A2,
            se=StdErr, P=Pvalue)]
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

info <- data.table(samples=1030836, cases=60620, controls=970216)
fwrite(info, sep="\t", quote=FALSE, file=sprintf("%s/%s/info.txt", out_dir, name))



