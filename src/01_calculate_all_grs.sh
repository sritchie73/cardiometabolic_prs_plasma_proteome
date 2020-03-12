#!/bin/bash

memory=25000 # set here because it also needs to be passed to plink

if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch." 
   exit 1
fi

# Before processing the genotype data, check whether there's any GRSs to calculate
nothing_to_do=1
for GRS_dir in data/GRS_resources/*; do
  GRS=$(basename $GRS_dir)
  if [ ! -d "analyses/GRS_profiles/$GRS" ]; then
    nothing_to_do=0
  fi
done

if [ "$nothing_to_do" -eq 1 ]; then
  echo "No new GRSs to calculate, aborting."
  exit 1
fi


# allows user to specify a job to wait for completion before running any of these scripts
if [ ! -z "$1" ]; then
  previous_job=$1
else 
  previous_job=1 # run first job immediately
fi

# First, we'll apply all the necessary filtering and generate plink-friendly 
# genotype files. This makes the bulk of the runtime when calculating any given
# GRS, so doing it upfront reduces runtime for any individual GRS.
wait_type="afterany" # for job after this one
if [ ! -d "analyses/processed_genotypes" ]; then
  mkdir -p analyses/processed_genotypes
  mkdir -p logs/genotype_qc
  qc_job=$(sbatch --dependency afterany:$previous_job \
                  --parsable \
									--account INOUYE-SL3-CPU \
                  --job-name "GRS qc INTERVAL" \
					   		  --array 1-22 \
                  --time 3:0:0 \
                  --mem $memory \
                  --output logs/genotype_qc/%A_genotype_qc_chr%a.stdout \
                  --error logs/genotype_qc/%A_genotype_qc_chr%a.stderr \
                  --partition skylake-himem \
                  src/01_job_scripts/01_genotype_qc.sh $memory)
  previous_job=$qc_job
  wait_type="afterok"
fi

# Iterates through GRS directories and checks if a corresponding directory 
# exists with the calculated profiles. If not, submits a batch script to
# calculate the profile. To avoid spawning too many jobs, each individual
# GRS profile calculation will wait for all others to finish.
#
# Calculating a GRS/PRS's levels in INTERVAL requires three consecutive jobs:
#
#  (1) To create a file of GRS weights were the effect allele has been
#      flipped.
#  (2) For each chromosome, sum the effect allele count * summary statistic
#      from GWAS for the trait of interest. This is done twice: once per
#      "strand".
#  (3) Combines all the per-chromosome files by summing the per-chromosome-per-strand
#      sums calculated in step  to obtain a total GRS/PRS level in each individual.
#      These GRS/PRS levels are then standardised (mean = 0, sd = 1) across INTERVAL.
#
# The reason for the strand flipping is that plink --score matches variants based on
# rsID and effect allele. If the GRS has been trained in a dataset that has been 
# genotyped or imputed to the opposite strand to INTERVAL, then those variants will 
# not be matched. plink will not try strand flipping if it finds a variant (by rsid)
# but with mismatching alleles (I have filed an issue requesting this feature:
# https://github.com/chrchang/plink-ng/issues/81). To maximise GRS coverage, we 
# therefore calculate the sum of effect alleles * summary stats twice: once matching
# the variants + effect alleles as-is, then again matching any missing variants with the 
# flipped effect allele to capture any cases where that variant is on the opposite strand
# in the INTERVAL genotype data. 
for GRS_dir in data/GRS_resources/*/; do
  GRS=$(basename $GRS_dir)
  if [ ! -d "analyses/GRS_profiles/$GRS" ]; then
    mkdir -p logs/GRS_profiles/$GRS
    mkdir -p analyses/GRS_profiles/$GRS
    # Create strand flipped GRS weights and recode variants to common
    # chr:pos:A1:A2 format (where A1:A2 alphabetical order)
    prep_job=$(sbatch --dependency $wait_type:$previous_job \
                      --kill-on-invalid-dep yes \
                      --parsable \
											--account INOUYE-SL3-CPU \
                      --job-name "$GRS preprocess" \
                      --time 1:0:0 \
                      --output logs/GRS_profiles/$GRS/%j_grs_preprocess.stdout \
                      --error logs/GRS_profiles/$GRS/%j_grs_preprocess.stderr \
                      --partition skylake \
											--cpus-per-task 1 \
                      src/01_job_scripts/02_grs_preprocess.sh $GRS)
    # Calculate GRS profile in each chromosome
    grs_job=$(sbatch --dependency afterok:$prep_job \
                     --kill-on-invalid-dep yes \
	 	  						   --parsable \
										 --account INOUYE-SL3-CPU \
                     --job-name "$GRS GRS profile" \
                     --array 1-22 \
                     --mem $memory \
                     --time 3:0:0 \
                     --output logs/GRS_profiles/$GRS/%A_calculate_profile_chr%a.stdout \
                     --error logs/GRS_profiles/$GRS/%A_calculate_profile_chr%a.stderr \
                     --partition skylake-himem \
                     src/01_job_scripts/03_calculate_profile.sh $GRS $memory)
    # Combine scores across chromosomes
    cleanup_job=$(sbatch --dependency afterok:$grs_job \
                         --kill-on-invalid-dep yes \
                         --parsable \
										     --account INOUYE-SL3-CPU \
												 --job-name "$GRS combine chr profile files" \
												 --time 1:0:0 \
												 --output logs/GRS_profiles/$GRS/%j_combine_chr_files.stdout \
												 --error logs/GRS_profiles/$GRS/%j_combine_chr_files.stderr \
												 --partition skylake \
												 --cpus-per-task 1 \
												 --ntasks 1 \
												 src/01_job_scripts/04_chr_profile_combiner.sh $GRS)
    echo "Submitted job to calculate $GRS GRS profile, assigned job ids: $prep_job,$grs_job,$cleanup_job"
    #previous_job=$grs_job # wait until computation intense job finishes before spawning the next
    wait_type="afterany"
  fi
done

# To save cluster space, remove the processed genotypes directory
#sbatch --dependency=afterany:$cleanup_job \
#       --time=1:0:0 \
#       --job-name=cleanup \
#       --wrap "rm -rf analyses/processed_genotypes"
