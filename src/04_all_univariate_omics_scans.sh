#!/bin/bash


# The purposes of this script is to submit a SLURM job that runs
# univariate association scans between each GRS and each omic platform.
# This pipeline is designed to be scalable to addition of GRSs and
# omic data as they become available, so first this script determines
# how many possible jobs there are to create one big array job.

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

# Create two arrays, one containing the GRS profile for each job,
# the other containing the omic profile for each job, as well as a
# counter tracking the total number of jobs to declare in the SLURM
# job array.
counter=0
declare -a "GRS_array"
declare -a "platform_array"

for GRS_dir in analyses/GRS_profiles/*/; do
  GRS=$(basename $GRS_dir)
  for omic_dir in analyses/processed_traits/*/; do
  omic_platform=$(basename $omic_dir)
    if [ ! -d "analyses/univariate_associations/$GRS/$omic_platform" ]; then
      GRS_array[$counter]=$GRS
      platform_array[$counter]=$omic_platform
      let "counter++"
    fi
  done
done

if [ "$counter" -eq "0" ]; then
  echo "Nothing to do, all association scans already run."
  exit 0
fi

mkdir -p "analyses/univariate_associations"
mkdir -p "logs/univariate_scans"
sbatch --dependency afterany:$previous_job \
       --job-name "GRS association scan" \
       --array 1-$counter \
       --time 1:0:0 \
       --output logs/univariate_scans/%A_%a_univariate_scan.stdout \
       --error logs/univariate_scans/%A_%a_univariate_scan.stderr \
       --partition cardio \
			 --account CARDIO-SL0-CPU \
       src/04_job_scripts/univariate_scan_dispatch.sh $counter "${GRS_array[@]}" "${platform_array[@]}"

