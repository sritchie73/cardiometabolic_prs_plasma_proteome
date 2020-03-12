#!/bin/bash

arg_array=($@)
total_jobs=${arg_array[0]}

GRS=${arg_array[$SLURM_ARRAY_TASK_ID]}
omic_platform=${arg_array[$SLURM_ARRAY_TASK_ID+$total_jobs]}

mkdir -p "analyses/univariate_associations/$GRS/$omic_platform"
echo $GRS $omic_platform | tee /dev/stderr
Rscript src/04_job_scripts/helpers/univariate_association_scan.R $GRS $omic_platform
