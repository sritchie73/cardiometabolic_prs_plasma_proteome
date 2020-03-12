#!/bin/bash

view_file=$1
grs=$(tail -n +2 $view_file | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f 1)

grs_script="src/07_job_scripts/01_helpers/${grs}.R"

if [ ! -f $grs_script ]; then
  echo "Script for extracting GWAS summary statistics for the $grs GRS has not been implemented!"
  exit 1
fi

ls analyses/mendelian_randomisation/$grs/gwas_summary_stats/*/info.txt &> /dev/null
if [ $? -eq 0 ]; then
  echo "already run, nothing to do"
  exit 0
fi

Rscript $grs_script

