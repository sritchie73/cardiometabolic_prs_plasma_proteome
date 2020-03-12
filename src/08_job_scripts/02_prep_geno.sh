#!/bin/bash

out_dir=analyses/grs_pqtl_removed/

for chr in $(seq 1 22); do 
  cp analyses/processed_genotypes/impute_chr$chr\_interval_filtered.bim $out_dir/chr$chr\.bim
  ln -s $PWD/analyses/processed_genotypes/impute_chr$chr\_interval_filtered.bed $out_dir/chr$chr\.bed
  ln -s $PWD/analyses/processed_genotypes/impute_chr$chr\_interval_filtered.fam $out_dir/chr$chr\.fam
  Rscript src/08_job_scripts/02_helpers/handle_duplicates.R $out_dir/chr$chr\.bim
done
