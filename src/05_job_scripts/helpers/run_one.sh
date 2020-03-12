#!/bin/bash

arg_array=($@)
total_jobs=${arg_array[0]}

GRS=${arg_array[$SLURM_ARRAY_TASK_ID]}
omic_platform=${arg_array[$SLURM_ARRAY_TASK_ID+$total_jobs]}
analysis_script=${arg_array[$total_jobs+$total_jobs+1]}
out_dir=${arg_array[$total_jobs+$total_jobs+2]}

echo $GRS $omic_platform | tee /dev/stderr
Rscript $analysis_script $GRS $omic_platform $out_dir
