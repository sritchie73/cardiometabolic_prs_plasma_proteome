#!/bin/bash

if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch."
   exit 1
fi

# Two optional arguments: view file of GRSs to restrict analysis to, and
# job id to wait for in the queue
if [ -z "$1" ]; then
  echo "You must provide a view file!"
  exit 1
else
  view_file=$1
fi

if [ ! -z "$2" ]; then
  previous_job=$2
else
  if [ ! -f $1 ]; then
		echo "You must provide a view file!"
		exit 1
  else
    previous_job=1
  fi
fi

# arguments wrong way round, swap
if [ -f $previous_job ]; then
  old_1=$view_file
  old_2=$previous_job
  previous_job=$old_1
  view_file=$old_2
fi

echo $previous_job
echo $view_file

mkdir -p analyses/mendelian_randomisation/
mkdir -p logs/mendelian_randomisation/

# Process disease GWAS summary statistics for each GRS:
n_grss=$(tail -n +2 $view_file | wc -l | cut -f 1 -d " ")
## ss_job=$(sbatch --dependency afterany:$previous_job \
##                 --parsable \
##                 --job-name "GWAS SS" \
##                 --time 1:0:0 \
##                 --array 1-$n_grss \
##                 --mem 16000 \
##                 --output logs/mendelian_randomisation/process_ss_%A_%a.out \
##                 --error logs/mendelian_randomisation/process_ss_%A_%a.err \
## 								--account INOUYE-SL2-CPU \
##                 --partition skylake \
##                 src/07_job_scripts/01_process_summary_stats.sh $view_file)
ss_job=$previous_job


# Call cis-pQTLs from full summary stats. This outputs two sets of results:
#
# (1) all cis-SNPs with P < 1e-4 for each protein
# (2) cis-SNPs where locally corrected P-value (Bonferroni) < threshold 
#     determined by FDR correction of sentinel SNP across all proteins
#
# The cis window is defined as 1MB upstream/downstream of the protein's gene's
# transcription start site. Protein complexes (i.e. multiple genes) are excluded. 
# 
# The latter set of results is based on recommendations for eQTL analysis from
#
# > Huang QQ et al.  Power, false discovery rate and Winner’s Curse in eQTL studies. 
#   Nucleic Acids Res. (2018). doi:10.1093/nar/gky780
#
# However, it is not clear whether this hierarchical approach will still control for FDR
# when looking at just a subset of genes (i.e. proteins on the somalogic panel). Bonferroni
# correction within each cis-window was also shown to be too conservative, and this may be
# exarcerbated by the increased density of the UK10K imputation panel. However, it is too
# computationally expensive to estimate the number of independent tests in each cis-window
# as that requires calculation of pairwise LD. 
#
# Note that the hierarchical correction is always more conservative than the 1e-4 threshold.
#
if [ ! -f "analyses/mendelian_randomisation/pqtls/cis_P_1e4.tsv" ]; then
  mkdir -p analyses/mendelian_randomisation/pqtls/ 
	cis_job=$(sbatch --dependency=afterany:$ss_job \
									 --parsable \
									 --job-name "cis-pQTLs" \
									 --time 2:0:0 \
									 --mem 16384 \
									 --output logs/mendelian_randomisation/call_cis_pqtls_%j.out \
									 --error logs/mendelian_randomisation/call_cis_pqtls_%j.err \
									 --account INOUYE-SL2-CPU \
									 --partition skylake \
									 --wrap "Rscript src/07_job_scripts/02_call_cis_pqtls.R")
else
  cis_job=$ss_job
fi

# For MR, we need to use cis-pQTLs, and preferably multiple (>3) *independent*
# cis-pQTLs so we can assess estimate robustness. To determine pQTL independence
# we will need a pairwise LD-map for all significant pQTLs that we can then later
# filter on the basis of the set of cis-pQTLs overlapping with each GWAS. We calculate
# the LD for all SNPs with P < 1e-4 since the hierarchical correction significant pQTLs
# are always a subset of these.
if [ ! -d "analyses/mendelian_randomisation/pqtls/SOMAMERs_cis_LD/" ]; then
  mkdir -p analyses/mendelian_randomisation/pqtls/SOMAMERs_cis_LD/
  nCores=16
	ld_job=$(sbatch --dependency afterany:$cis_job \
								  --parsable \
								  --job-name "cis-pQTL-LD" \
								  --time 12:0:0 \
								  --ntasks 1 \
								  --cpus-per-task $nCores \
								  --output logs/mendelian_randomisation/cis_pqtl_ld_%j.out \
								  --error logs/mendelian_randomisation/cis_pqtl_ld_%j.err \
									--account INOUYE-SL2-CPU \
								  --partition skylake,skylake-himem \
								  --wrap "Rscript src/07_job_scripts/03_cis_pqtl_ld.R $nCores")
else
  ld_job=$cis_job
fi

# For the published pQTLs (mix of cis and trans) we often have cases where the 
# pQTL is not in the GWAS, so we will need to use proxy variants in those cases 
# (closest SNP in high LD: r2 > 0.8 within 250kb).
if [ ! -d "analyses/mendelian_randomisation/pqtls/conditional_pQTLs_tag_LD" ]; then
  mkdir -p analyses/mendelian_randomisation/pqtls/conditional_pQTLs_tag_LD
  nCores=6
	ld_job2=$(sbatch --dependency afterany:$ld_job \
								   --parsable \
								   --job-name "trans-pQTL-LD" \
								   --time 10:0:0 \
								   --ntasks 1 \
								   --cpus-per-task $nCores \
								   --output logs/mendelian_randomisation/trans_pqtl_ld_%j.out \
								   --error logs/mendelian_randomisation/trans_pqtl_ld_%j.err \
									 --account INOUYE-SL2-CPU \
								   --partition skylake \
								   --wrap "Rscript src/07_job_scripts/04_published_pqtl_ld.R $nCores")
else
  ld_job2=$ld_job
fi

# For the PRSs of interest, map their pQTLs to their GWAS summary statistics
nPRS=$(sed 1d $view_file | wc -l)
map_job=$(sbatch --dependency afterany:$ld_job2 \
                 --parsable \
                 --job-name "Map summary stats" \
                 --time 24:0:0 \
                 --array 1-$nPRS \
                 --mem 64000 \
                 --output logs/mendelian_randomisation/map_summary_stats_%A_%a.out \
                 --error logs/mendelian_randomisation/map_summary_stats_%A_%a.err \
								 --account INOUYE-SL2-CPU \
								 --partition skylake,skylake-himem \
                 --wrap "Rscript src/07_job_scripts/05_map_summary_stats.R $view_file")

# Test for co-localisation of pQTLs with GWAS summary statistics
nPRS=$(sed 1d $view_file | wc -l)
coloc_job=$(sbatch --dependency=afterany:$map_job \
									 --parsable \
									 --job-name "GRS coloc" \
									 --time 5:0:0 \
									 --array 1-$nPRS \
									 --mem 16000 \
									 --output logs/mendelian_randomisation/colocalisation_%A_%a.out \
									 --error logs/mendelian_randomisation/colocalisation_%A_%a.err \
									 --account INOUYE-SL2-CPU \
									 --partition skylake \
									 --wrap "Rscript src/07_job_scripts/06_colocalisation.R $view_file")

# Finally, run mendelian randomisation analysis
nPRS=$(sed 1d $view_file | wc -l)
mr_job=$(sbatch --dependency afterany:$coloc_job \
								--parsable \
								--job-name "GRS MR" \
								--time 3:0:0 \
								--array 1-$nPRS \
								--mem 16000 \
								--output logs/mendelian_randomisation/mr_%A_%a.out \
								--error logs/mendelian_randomisation/mr_%A_%a.err \
							  --account INOUYE-SL2-CPU \
							  --partition skylake \
								--wrap "Rscript src/07_job_scripts/07_mendelian_randomisation.R $view_file")

