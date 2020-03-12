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

mkdir -p analyses/grs_pqtl_removed
mkdir -p logs/grs_pqtl_removed

# Check if another job is already running, if so abort, 
# otherwise get pairs to run
if [ -f analyses/grs_pqtl_removed/pairs_to_run.txt ]; then
  echo "Another job of this type is running, you will need to wait for that to finish."
  exit 1
else
	Rscript src/08_job_scripts/01_grs_apt_pairs.R $view_file
	pair_file=analyses/grs_pqtl_removed/pairs_to_run.txt
	nPairs=$(wc -l $pair_file | cut -f 1 -d " ")
fi

# Prepare genotype data
prep_job=$(sbatch --dependency afterany:$previous_job \
                  --parsable \
                  --job-name "GRS pQTL removal: process bim files" \
                  --time 1:0:0 \
                  --output logs/grs_pqtl_removed/prep_geno_%j.out \
                  --error logs/grs_pqtl_removed/prep_geno_%j.err \
									--account INOUYE-SL3-CPU \
									--partition skylake \
                  src/08_job_scripts/02_prep_geno.sh)
 

# Extract pQTL variants to remove
extract_job=$(sbatch --dependency afterok:$prep_job \
									   --kill-on-invalid-dep yes \
									   --parsable \
                     --job-name "pQTL tags" \
									   --time 1:0:0 \
									   --output logs/grs_pqtl_removed/extract_pQTL_ld_%j.out \
									   --error logs/grs_pqtl_removed/extract_pQTL_ld_%j.err \
										 --account INOUYE-SL3-CPU \
										 --partition skylake-himem \
                     --wrap "Rscript src/08_job_scripts/03_get_pqtl_windows.R")

last_job=$extract_job
wait_type=afterok
final_wait_list=""
for pair in $(seq 1 $nPairs); do
  GRS=$(sed "${pair}q;d" $pair_file | cut -f 1)
  Aptamer=$(sed "${pair}q;d" $pair_file | cut -f 2)

  echo "$GRS $Aptamer"

	# Remove pQTL and tags from the genotype data prior to GRS re-calculation
	filter_job=$(sbatch --dependency $wait_type:$last_job \
											--parsable \
											--kill-on-invalid-dep yes \
											--array 1-22 \
											--job-name "Filter pQTL tags" \
											--time 1:0:0 \
											--output logs/grs_pqtl_removed/$GRS\_$Aptamer\_filter_pQTL_tags_%A_%a.out \
											--error logs/grs_pqtl_removed/$GRS\_$Aptamer\_filter_pQTL_tags_%A_%a.err \
											--account INOUYE-SL3-CPU \
											--partition skylake \
											src/08_job_scripts/04_filter_tags.sh $GRS $Aptamer) 

	# Create strand flipped GRS weights since plink can't handle this yet
	prep_job=$(sbatch --dependency afterok:$filter_job \
										--kill-on-invalid-dep yes \
										--parsable \
										--job-name "$GRS GRS strand flip" \
										--time 1:0:0 \
										--output logs/grs_pqtl_removed/$GRS\_$Aptamer\_%j_grs_strand_flip.stdout \
										--error logs/grs_pqtl_removed/$GRS\_$Aptamer\_%j_grs_strand_flip.stderr \
									  --account INOUYE-SL3-CPU \
										--partition skylake \
										src/08_job_scripts/05_grs_preprocess.sh $GRS $Aptamer)

	# Calculate GRS profile in each chromosome
	grs_job=$(sbatch --dependency afterok:$prep_job \
									 --kill-on-invalid-dep yes \
									 --parsable \
									 --job-name "$GRS GRS profile" \
									 --array 1-22 \
									 --time 3:0:0 \
                   --mem 25000 \
									 --output logs/grs_pqtl_removed/$GRS\_$Aptamer\_%A_calculate_profile_chr%a.stdout \
									 --error logs/grs_pqtl_removed/$GRS\_$Aptamer\_%A_calculate_profile_chr%a.stderr \
									 --account INOUYE-SL3-CPU \
									 --partition skylake \
									 src/08_job_scripts/06_calculate_profile.sh $GRS $Aptamer 8192)

	# Combine scores across chromosomes
	combine_job=$(sbatch --dependency afterok:$grs_job \
											 --kill-on-invalid-dep yes \
											 --parsable \
											 --job-name "$GRS combine chr profile files" \
											 --time 1:0:0 \
											 --output logs/grs_pqtl_removed/$GRS\_$Aptamer\_%j_combine_chr_files.stdout \
											 --error logs/grs_pqtl_removed/$GRS\_$Aptamer\_%j_combine_chr_files.stderr \
											 --account INOUYE-SL3-CPU \
											 --partition skylake \
											 src/08_job_scripts/07_chr_profile_combiner.sh $GRS $Aptamer)

  # No need to throttle job throughput on CSD3
  # last_job=$combine_job
  # wait_type=afterany

  # But we instead need to collate a big list of jobs to wait for
  if [ $pair -eq 1 ]; then
    final_wait_list=$combine_job
  else
    final_wait_list="$final_wait_list,$combine_job"
  fi
done

# Re-test association
test_job=$(sbatch --dependency afterok:$final_wait_list \
									--kill-on-invalid-dep yes \
									--parsable \
                  --job-name "GRS test no pqtl" \
                  --time 1:0:0 \
                  --output logs/grs_pqtl_removed/%j_pair_retest.stdout \
                  --error logs/grs_pqtl_removed/%j_pair_retest.stderr \
									--account INOUYE-SL3-CPU \
									--partition skylake \
                  --wrap "Rscript src/08_job_scripts/08_retest.R $pair_file")

# Final cleanup
combine_job=$(sbatch --dependency afterok:$test_job \
										 --kill-on-invalid-dep yes \
										 --parsable \
										 --job-name "pQTL GRS exclusion cleanup" \
										 --time 0:10:0 \
										 --output logs/grs_pqtl_removed/final_cleanup_%j.out \
										 --error logs/grs_pqtl_removed/final_cleanup_%j.err \
										 --account INOUYE-SL3-CPU \
										 --partition skylake \
										 src/08_job_scripts/09_cleanup.sh)
