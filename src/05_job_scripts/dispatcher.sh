#!/bin/bash

# Set by top level script
if [ -z "$master_script" ]; then
  echo "This script should not be executed directly."
  exit 1
fi

view_file=$1
previous_job=$2
out_name=$3
analysis_script=$4

out_dir=analyses/sensitivity_associations/$out_name
mkdir -p $out_dir

# Create two arrays, one containing the GRS profile for each job,
# the other containing the omic profile for each job, as well as a
# counter tracking the total number of jobs to declare in the SLURM
# job array.
counter=0
declare -a "GRS_array"
declare -a "platform_array"

if [[ "$view_file" == "all" ]]; then
  for GRS_dir in analyses/GRS_profiles/*/; do
    GRS=$(basename $GRS_dir)
    for omic_dir in analyses/processed_traits/*/; do
       omic_platform=$(basename $omic_dir)
       if [ ! -d "$out_dir/$GRS/$omic_platform" ]; then
         mkdir -p $out_dir/$GRS/$omic_platform
         GRS_array[$counter]=$GRS
         platform_array[$counter]=$omic_platform
         let "counter++"
       fi 
    done
  done
else 
  GRSs=$(cut -f 1 $view_file | tail -n +2)
  for GRS in $GRSs; do 
    for omic_dir in analyses/processed_traits/*/; do
       omic_platform=$(basename $omic_dir)
       if [ ! -d "$out_dir/$GRS/$omic_platform" ]; then
         mkdir -p $out_dir/$GRS/$omic_platform
         GRS_array[$counter]=$GRS
         platform_array[$counter]=$omic_platform
         let "counter++"
       fi
    done
  done
fi

if [ "$counter" -eq "0" ]; then
  echo "Nothing to do, all association scans already run."
  exit 0
fi

mkdir -p "logs/sensitivity_associations/"
sbatch --dependency=afterany:$previous_job \
       --job-name="$out_name association scan" \
       --array=1-$counter \
       --mem-per-cpu=2048 \
       --time=0:10:0 \
       --output=logs/sensitivity_associations/%A_%a_$out_name\_scan.stdout \
       --error=logs/sensitivity_associations/%A_%a_$out_name\_scan.stderr \
       --account INOUYE-SL3-CPU \
       --partition skylake \
       src/05_job_scripts/helpers/run_one.sh $counter "${GRS_array[@]}" "${platform_array[@]}" $analysis_script $out_dir
