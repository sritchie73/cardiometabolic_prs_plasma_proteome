#!/bin/bash

if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch." 
   exit 1
fi

# allows user to specify a job to wait for completion before running any of these scripts
if [ ! -z "$1" ]; then
  previous_job=$1
else
  previous_job=1 # run first job immediately
fi

# Curate the omic platform files into a common format for polygenic
# association scans. Each should have a folder containing 'traits.tsv'
# with three columns: IID, variable, and value. Where IID is the genetic
# identifier (samples are subset to those with genetic information), 
# variable is a measurement identifier on the platform (e.g. a single protein), 
# and value is the value of that measurement in that sample.
#
# Additional files containing per-variable and per-sample information are
# also created, as well as platform-specific covariate files.
#
# A phenotypes.tsv file is also created in analyses/processed_traits/
# containing sample information, e.g. age, sex, etc. and clinical biomarkers
# not used in the association scan.
prep_job=$(sbatch --dependency afterany:$previous_job \
             		  --job-name "data prep" \
                  --time 3:0:0 \
                  --parsable \
									--account INOUYE-SL3-CPU \
                  --output logs/prep_omics_dirs_%j.stdout \
                  --error logs/prep_omics_dirs_%j.stderr \
                  --partition skylake \
									--ntasks 1 \
									--cpus-per-task 10 \
                  src/03_job_scripts/01_process_traits.sh)

# Process QTL files so that we have a list of variables and variants in the omic dir
qtl_job=$(sbatch --dependency aftercorr:$prep_job \
                 --kill-on-invalid-dep yes \
                 --parsable \
                 --job-name "process QTLs" \
                 --time 0:10:0 \
                 --output logs/process_qtls_%j.stdout \
                 --error logs/process_qtls_%j.stderr \
								 --account INOUYE-SL3-CPU \
								 --partition skylake \
								 --ntasks 1 \
								 --cpus-per-task 1 \
                 --wrap "Rscript src/03_job_scripts/03_process_qtls.sh")

# Extract genotype data so we can load it into R. 
geno_job=$(sbatch --dependency aftercorr:$qtl_job \
                  --kill-on-invalid-dep yes \
                  --parsable \
                  --array 1-22%6 \
                  --job-name "extract geno" \
                  --time 3:0:0 \
                  --output logs/extract_geno_%A_%a.stdout \
                  --error logs/extract_geno_%A_%a.stderr \
								  --account INOUYE-SL3-CPU \
								  --partition skylake-himem \
                  --mem 25000 \
                  src/03_job_scripts/04_extract_geno_prob.sh)

# Need to separately extract QTLs with only positional identifiers.
n_pos=$(grep ":" analyses/processed_traits/somalogic_proteins/qtl_variants.txt | wc -l | cut -f 1)
pos_job=$(sbatch  --dependency aftercorr:$geno_job \
                  --kill-on-invalid-dep yes \
                  --parsable \
                  --array 1-$n_pos%6 \
                  --job-name "extract positional geno" \
                  --time 3:0:0 \
                  --output logs/extract_geno_pos_%A_%a.stdout \
                  --error logs/extract_geno_pos_%A_%a.stderr \
								  --account INOUYE-SL3-CPU \
								  --partition skylake-himem \
                  --mem 25000 \
                  src/03_job_scripts/05_extract_geno_pos.sh)

# Combine the files
combine_job=$(sbatch --dependency aftercorr:$pos_job \
                     --kill-on-invalid-dep yes \
                     --parsable \
                     --job-name "combine geno" \
                     --time 1:0:0 \
                     --output logs/combine_geno_%j.stdout \
                     --error logs/combine_geno_%j.stderr \
									   --partition skylake \
									   --ntasks 1 \
									   --cpus-per-task 1 \
                     src/03_job_scripts/06_combine_geno.sh)
