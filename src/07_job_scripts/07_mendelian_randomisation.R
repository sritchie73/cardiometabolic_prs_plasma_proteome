library(data.table)
library(foreach)
library(MendelianRandomization)
library(scales)
library(ggplot2)
source("src/07_job_scripts/07_helpers/mr_functions.R")

# Determine GRS we're working with
if (!exists("grs")) {
	view_file <- commandArgs(trailingOnly=TRUE)[1]
	view <- fread(view_file)
	task <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
	grs <- view[task, GRS_name]

	# Paths:
	out_dir <- sprintf("analyses/mendelian_randomisation/%s", grs)

	# Check if there is anything to run:
	if (!file.exists(sprintf("%s/instruments.tsv", out_dir))) {
		quit(save="no", status=0)
	}
	if (dir.exists(sprintf("%s/dose_response_curves", out_dir))) {
		quit(save="no", status=0)
	}
} else {
	# Paths:
	out_dir <- sprintf("analyses/mendelian_randomisation/%s", grs)
}

# Load instruments:
instruments <- fread(sprintf("%s/instruments.tsv", out_dir))

# Drop instruments we couldn't map to GWAS summary stats:
instruments <- instruments[!is.na(effect.gwas)]

# Filter to protein level
instruments <- unique(instruments, by=c("GWAS", "PRS", "Target", "UniProt", "Gene", "chr", "pos", "EA", "OA"))

# Restrict to proteins with at least three IVs, and at least one cis-
iv_stats <- instruments[,.(N=.N, cis=any(type == "cis")), by=.(PRS, GWAS, Target, UniProt, Gene)]
pass <- iv_stats[N >= 3 & (cis), .(PRS, GWAS, Target, UniProt, Gene)]
instruments <- instruments[pass, on=.(PRS, GWAS, Target, UniProt, Gene)]

if (nrow(instruments) == 0) {
	if (!exists("grs")) {
    quit(save="no", status=0)
  } else {
    stop(grs)
  }
}

# Run MR analysis for each protein and target GWAS:
mr_dt <- foreach(gwas = unique(instruments$GWAS), .combine=rbind, .errorhandling="remove") %do% {
  prots <- unique(instruments[,.(PRS, GWAS, Target, UniProt, Gene)])[GWAS == gwas]
  foreach(protIdx = prots[,.I], .combine=rbind, .errorhandling="remove") %do% {
    prot <- prots[protIdx]

    # Filter tables to this exposure:
    sub_inst <- instruments[prot, on = .(PRS, GWAS, Target, UniProt, Gene)]

    # Construct MR input object:
		mri <- mr_input(
			bx = sub_inst$Prot.Effect.pQTL,
			bxse = sub_inst$Prot.SE.pQTL,
			by = sub_inst$effect.gwas,
			byse = sub_inst$se.gwas,
			snps = sub_inst$ALT_ID,
			effect_allele = sub_inst$EA,
			other_allele = sub_inst$OA,
			eaf = sub_inst$EAF.pqtl,
			exposure = prot$Gene,
			outcome = gwas)
    
    # Run all available MR methods:
    mr_results <- mr_tryall(mri)

    # Get GWAS type to appropriately label Y-axes
    gwas_info <- fread(sprintf("%s/gwas_summary_stats/%s/info.txt", out_dir, gwas))
    if ("cases" %in% names(gwas_info)) {
      ytitle <- sprintf("Log Odds of %s", gwas)
    } else {
      ytitle <- sprintf("SD effect on %s", gwas)
    }

    # Generate dose response curve
    plot_dir <- sprintf("%s/dose_response_curves/%s", out_dir, gwas)
    dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
    g <- gg_dose_response(sub_inst, mr_results, label=FALSE) +
         ylab(ytitle) +
         xlab(sprintf("SD effect on %s", prot$Gene))
    ggsave(g, file=sprintf("%s/%s.pdf", plot_dir, prot$Gene), width=8, height=7.2, useDingbats=FALSE)

    # Add annotations
    mr_results <- cbind(prot, mr_results)

    # Return results
    return(mr_results)
  }
}
fwrite(mr_dt, sep="\t", quote=FALSE, file=sprintf("%s/mr_results.tsv", out_dir))


