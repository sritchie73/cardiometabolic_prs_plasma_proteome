#!/bin/bash

for pt_dir in analyses/processed_traits/*; do
  if [ -f $pt_dir/qtl_list.txt ] && [ ! -f $omic_dir/qtl_geno_prob.tsv ]; then
    Rscript src/03_job_scripts/06_helpers/combine_geno.R $pt_dir
    rm $pt_dir/*.raw
    rm $pt_dir/*.log
    rm $pt_dir/*.nosex
    rm $pt_dir/*.bed
    rm $pt_dir/*.bim
    rm $pt_dir/*.fam
    rm $pt_dir/*.map
    rm $pt_dir/*.ped
  fi
done
