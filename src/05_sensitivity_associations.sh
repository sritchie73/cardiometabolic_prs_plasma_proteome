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


# Get the maximum job id of all jobs running for a user -
# allows us to dispatch multiple job scripts
last_active_job () {
  # get the largest job id for active jobs in the queue for the current user
  max_job_id=$(squeue -u $(whoami) -h -o %F | sort -nr | head -n 1)

  # If no jobs in queue, return 1 (run immediately) 
  if [ -z $max_job_id ]; then
    max_job_id=1
  fi

  # If there are jobs in the queue, but its the same max job ID 
  # as before running any dispatch scripts (e.g. script exited
  # with nothing to run), return the original $previous_job
  # specified by the user:
  if [ ! -z $current_max_job ]; then
     if [ $current_max_job -eq $max_job_id ]; then
       max_job_id=$previous_job
     fi
  fi

  echo $max_job_id
}

# variable to put in the child process environment so we
# can check the next job scripts are being executed directly
export master_script=1

current_max_job=$(last_active_job)

# Run association analyses adjusted for BMI
bash src/05_job_scripts/dispatcher.sh \
   $view_file \
   $previous_job \
   BMI_adjusted \
   src/05_job_scripts/analysis_scripts/univariate_adjusted_for_BMI.R

# Run association analyses adjusting for QTLs
bash src/05_job_scripts/dispatcher.sh \
	$view_file \
  $(last_active_job) \
  qtl_prob_dosage_adjusted \
  src/05_job_scripts/analysis_scripts/univariate_adjusted_for_qtls_prob_dosage.R

# Adjust associations for circadian and season effects
bash src/05_job_scripts/dispatcher.sh \
  $view_file \
  $(last_active_job) \
  circadian_adjusted \
  src/05_job_scripts/analysis_scripts/univariate_adjusted_for_circadian.R

bash src/05_job_scripts/dispatcher.sh \
  $view_file \
  $(last_active_job) \
  season_adjusted \
  src/05_job_scripts/analysis_scripts/univariate_adjusted_for_season.R

# Run view-specific adjustments
if [ $view_file == "views/lung_function.txt" ]; then 
	bash src/05_job_scripts/dispatcher.sh \
		$view_file \
		$(last_active_job) \
		height_adjusted \
		src/05_job_scripts/analysis_scripts/univariate_adjusted_for_height.R
fi





