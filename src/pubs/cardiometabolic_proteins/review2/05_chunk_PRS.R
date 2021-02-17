suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(R.utils))

out_dir="analyses/pub/cardiometabolic_proteins/review2"

# Setup local environment
chrIdx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ncores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
if(is.na(ncores)) ncores <- 1
setDTthreads(ncores)

# fread can fall over when concurrently attempting to read the same score file from 22 processes.
# This function attempts to make this "task safe" by allowing only one call to fread to read from
# a given file at a time
sfread <- function(file, ...) {
  fname <- basename(dirname(file))
  lockfile <- sprintf("%s/00FREADLOCK-%s-grs-weights", out_dir, fname)
  while(file.exists(lockfile)) {
    Sys.sleep(5)
  }
  cat(chrIdx, "\n", file=lockfile)
  on.exit({ system(sprintf("rm -f %s", lockfile), wait=TRUE) })
  fread(file, ...)
}

# Load in three scores and filter to this chromosome
scores <- rbind(idcol="PRS",
  AF_PRS=sfread("data/GRS_resources/Afib_2018/grs_weights.txt.gz"),
  CAD_PRS=sfread("data/GRS_resources/CAD_metaGRS/grs_weights.txt.gz"),
  CKD_PRS=sfread("data/GRS_resources/CKD_2019/grs_weights.txt.gz"),
  IS_PRS=sfread("data/GRS_resources/Stroke_metaGRS/grs_weights.txt.gz"),
  T2D_PRS=sfread("data/GRS_resources/T2D_2018/grs_weights.txt.gz")
)
scores <- scores[chr == chrIdx]

# Drop strand ambiguous variants:
flip_strand <- function(x) {
  # Swap each letter for a dummy, we need this intermediate
  # step so we can distinguish between alleles when swapping.
  # E.g if we did A -> T then T -> A we'd end up with all A's
  # and no T's. instead we do A -> V -> T and T -> X -> A.
  x <- gsub("A", "V", x)
  x <- gsub("T", "X", x)
  x <- gsub("C", "Y", x)
  x <- gsub("G", "Z", x)
  x <- gsub("V", "T", x)
  x <- gsub("X", "A", x)
  x <- gsub("Y", "G", x)
  x <- gsub("Z", "C", x)
  return(x)
}
scores <- scores[effect_allele != flip_strand(other_allele)]

# Load in BIM file and orient effect alleles to the correct strand.
bim <- fread(sprintf("analyses/processed_genotypes/impute_chr%s_interval_filtered.bim", chrIdx), header=FALSE)
setnames(bim, c("chr", "rsid", "cm", "pos", "A1", "A2"))

scores[, match := FALSE]
scores[bim, on = .(pos, effect_allele=A1, other_allele=A2), match := TRUE]
scores[bim, on = .(pos, effect_allele=A2, other_allele=A1), match := TRUE]

scores[bim, on = .(pos), c("effect_allele", "other_allele") := 
       .(ifelse(match, effect_allele, flip_strand(effect_allele)),
         ifelse(match, other_allele, flip_strand(other_allele)))]
scores[bim, on = .(pos, effect_allele=A1, other_allele=A2), match := TRUE]
scores[bim, on = .(pos, effect_allele=A2, other_allele=A1), match := TRUE]

# Drop variants without any match:
scores <- scores[(match)]

# Give each variant a unique identifier
bim[, rsid := paste(chr, pos, A1, A2, sep=":")]
scores[bim, on = .(pos, effect_allele=A1, other_allele=A2), rsid := i.rsid]
scores[bim, on = .(pos, effect_allele=A2, other_allele=A1), rsid := i.rsid]

# Drop any multi-allelic score variants
mult <- scores[,.N,by=.(PRS, pos)][N > 1]
scores <- scores[!mult, on = .(PRS, pos)]

# Orient variants to the same effect allele (plink warns or crashes otherwise)
scores[bim, on = .(pos, effect_allele=A2, other_allele=A1), 
       c("effect_allele", "other_allele", "weight") := 
       .(other_allele, effect_allele, -weight)]

# Partition variants into independend LD blocks estimated by Berisa et al. 2016
# PMID: 26395773
# Data: https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-all.bed
ld_blocks <- fread("data/Berisa_etal_2016/EUR_1000G_ind_ld_blocks.bed")
ld_blocks[, block_num := .I]
ld_blocks[, chr := as.integer(gsub("chr", "", chr))]

scores[ld_blocks, on = .(chr, pos >= start, pos < stop), block_num := i.block_num]
scores[pos == ld_blocks[chr == chrIdx, max(stop)], block_num := ld_blocks[chr == chrIdx, max(block_num)]]

# spread to wide format:
scores[, block_num := paste0(PRS, "_", block_num)]
scores <- dcast(scores, rsid + effect_allele ~ block_num, value.var="weight", fill=0)

# Write out score file and variant information with new IDs and symlink to
# remaining genotype data
fwrite(scores, sep="\t", quote=FALSE, file=sprintf("%s/scores_chr%s.txt", out_dir, chrIdx))
fwrite(scores[,.(rsid)], quote=FALSE, col.names=FALSE, file=sprintf("%s/variants_chr%s.txt", out_dir, chrIdx))
fwrite(bim, sep="\t", quote=FALSE, col.names=FALSE, sprintf("%s/chr%s.bim", out_dir, chrIdx))
system(sprintf("ln -s $(realpath analyses/processed_genotypes/impute_chr%s_interval_filtered.bed) %s/chr%s.bed", chrIdx, out_dir, chrIdx), wait=TRUE)
system(sprintf("ln -s $(realpath analyses/processed_genotypes/impute_chr%s_interval_filtered.fam) %s/chr%s.fam", chrIdx, out_dir, chrIdx), wait=TRUE)

# How many columns in the score file?
ncols <- ncol(scores)

# remove uncessary objects in memory to maximise memory available to plink2
rm(scores, bim, mult)
invisible(gc())

# Run plink2 --score
cmd="plink2"
cmd[2] <- sprintf("--bfile %s/chr%s", out_dir, chrIdx)
cmd[3] <- sprintf("--out %s/scores_chr%s", out_dir, chrIdx)
cmd[4] <- sprintf("--threads %s --memory 23488 --silent", ncores)
cmd[5] <- sprintf("--extract %s/variants_chr%s.txt", out_dir, chrIdx)
cmd[6] <- sprintf("--score %s/scores_chr%s.txt", out_dir, chrIdx)
cmd[7] <- sprintf("'header-read' 'ignore-dup-ids' 'cols=scoresums'")
cmd[8] <- sprintf("--score-col-nums 3-%s", ncols)
system(paste(na.omit(cmd), collapse=" "), wait=TRUE)

# Create file to let others know this chromosome has finished processing
system(sprintf("touch %s/00FIN-chr-%s", out_dir, chrIdx), wait=TRUE)

# Last task to finish collates the results.
finfiles <- list.files(path=out_dir, pattern="00FIN-chr-*")
if (length(finfiles) <  22) {
  quit(save="no") # other tasks still running, finish.
}

lockfile <- sprintf("%s/00LOCK-collate", out_dir)
if (!file.exists(lockfile)) {
  # If two or more tasks reach this point at the same time, each one sleeps for a different
  # amount of time then checks for lock file again to prevent multiple tasks doing the collation
  # step and messing up the output.
  Sys.sleep(chrIdx)
}

if (file.exists(lockfile)) {
  quit(save="no") # Another task that finished at the same time and has taken over the collation process
}

# Create lock file and enter the task ID for debugging purposes.
cat(chrIdx, "\n", file=lockfile)

# remove no longer needed 00FIN files
for (chr in 1:22) {
  system(sprintf("rm -f %s/00FIN-chr-%s", out_dir, chr), wait=TRUE)
}

# Collate plink log files
logfile <- sprintf("%s/collated_plink_logs.txt", out_dir)
system(sprintf("touch %s", logfile), wait=TRUE)
for (chr in 1:22) {
  system(sprintf("cat %s/scores_chr%s.log >> %s", out_dir, chr, logfile), wait=TRUE)
  system(sprintf("rm -f %s/scores_chr%s.log", out_dir, chr), wait=TRUE)
}

# Collate sscore files (each chromosome has different "scores")
inner_join <- function(a, b) { a[b, on = .(IID), nomatch=0] }
sscores <- foreach(chr = 1:22, .combine=inner_join) %do% {
  dt <- fread(sprintf("%s/scores_chr%s.sscore", out_dir, chr))
  setnames(dt, "#IID", "IID")
  system(sprintf("rm -f %s/scores_chr%s.sscore", out_dir, chr), wait=TRUE)
  system(sprintf("rm -f %s/scores_chr%s.txt", out_dir, chr), wait=TRUE)
  system(sprintf("rm -f %s/variants_chr%s.txt", out_dir, chr), wait=TRUE)
  return(dt)
}
names(sscores) <- gsub("_SUM", "", names(sscores))
sscores <- melt(sscores, id.vars="IID", variable.name="score_chunk", value.name="score_sum")
fwrite(sscores, sep="\t", quote=FALSE, file=sprintf("%s/chunked_scores.sscore.gz", out_dir), compress="gzip")

# Clean up remaining files
for (chr in 1:22) {
  system(sprintf("rm -f %s/chr%s.bim", out_dir, chr), wait=TRUE)
  system(sprintf("rm -f %s/chr%s.bed", out_dir, chr), wait=TRUE)
  system(sprintf("rm -f %s/chr%s.fam", out_dir, chr), wait=TRUE)
}
system(sprintf("rm -f %s/00LOCK-collate", out_dir), wait=TRUE)
system(sprintf("rm -rf %s/slurm_logs/", out_dir), wait=TRUE)

