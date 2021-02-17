library(data.table)
library(foreach)
library(doMC)

# Set up parallel environment
ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
ncores <- as.integer(ncores)
if(is.na(ncores)) ncores <- 1
registerDoMC(ncores)
setDTthreads(ncores)

out_dir <- "analyses/mendelian_randomisation/pqtls/"

# Load in protein information and split out protein complexes and
# proteins encoded by multiple genes
prots <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
locmap <- prots[,.(
  Gene=strsplit(Gene.Name, "\\|")[[1]],
  chr=strsplit(chr, "\\|")[[1]],
  start=strsplit(start, "\\|")[[1]]),
  by=.(SOMAMER_ID, SeqId)]
locmap <- locmap[,.(start=strsplit(start, ";")[[1]]), by=.(SOMAMER_ID, SeqId, Gene, chr)]
locmap <- locmap[chr %in% 1:22]

cis_stats <- foreach(somamer=unique(locmap$SOMAMER_ID), .combine=rbind, .errorhandling="remove") %dopar% {
  locs <- locmap[SOMAMER_ID == somamer]
  chr_stats <- foreach(lIdx = seq_len(locs[,.N]), .combine=rbind, .errorhandling="remove") %do% {
    chr <- locs[lIdx, chr]
    tss <- as.numeric(locs[lIdx, start])
    loc_stats <- fread(sprintf("data/full_pQTL_summary_stats/%s/%s_chrom_%s_meta_1.tbl.gz", somamer, somamer, chr))
    loc_stats <- loc_stats[position > tss - 1e6 & position < tss + 1e6]
    loc_stats
  }
  chr_stats <- unique(chr_stats)
  chr_stats[, .(SOMAMER_ID=somamer, chr=chromosome, pos=position, 
                EA=toupper(Allele1), OA=toupper(Allele2), 
                effect=Effect, se=StdErr, P=10^(`log(P)`))]
}

fwrite(cis_stats[P < 1e-4], sep="\t", quote=FALSE, file=sprintf("%s/cis_P_1e4.tsv", out_dir))

# Hierarchical correction:
cis_stats[, local_bonf := pmin(P*.N, 1), by=SOMAMER_ID]
global <- cis_stats[, .SD[which.min(P)], by=SOMAMER_ID]
global[, global_fdr := p.adjust(local_bonf, method="fdr")]
cis_stats[global, on = .(SOMAMER_ID), global_fdr := i.global_fdr]

# Determine local correction threshold corresponding to global fdr < 0.05
global <- global[order(global_fdr)]
sig_idx <- global[global_fdr < 0.05, .N]
threshold <- global[sig_idx:(sig_idx+1), mean(local_bonf)]

# Log p-value threshold for each aptamer
cis_thresh <- cis_stats[,.(threshold = threshold/.N),by=SOMAMER_ID]
fwrite(cis_thresh, sep="\t", quote=FALSE,
       file=sprintf("%s/aptamer_cis_hierarchical_pthresh.tsv", out_dir))

# Filter
fwrite(cis_stats[local_bonf < threshold], sep="\t", quote=FALSE, 
       file=sprintf("%s/cis_hierarchical_correction.tsv", out_dir))



