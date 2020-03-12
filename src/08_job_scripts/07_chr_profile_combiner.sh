#!/bin/bash

grs_name=$1
Aptamer=$2
out_folder=analyses/grs_pqtl_removed/$grs_name/$Aptamer
grs_folder=data/GRS_resources/$grs_name

# Combine all the score files across chromosomes
Rscript src/08_job_scripts/07_helpers/chr_profile_combiner.R $grs_name $Aptamer

# Combine all variant lists across chromosomes
for chr in $(seq 1 22); do
  cat $out_folder/flip_profile_chr$chr.sscore.vars >> $out_folder/strand_mismatch.sscore.vars
done
gzip $out_folder/strand_mismatch.sscore.vars

for chr in $(seq 1 22); do
   cat $out_folder/profile_chr$chr.sscore.vars >> $out_folder/profile.sscore.vars
   cat $out_folder/flip_profile_chr$chr.sscore.vars >> $out_folder/profile.sscore.vars
done
gzip $out_folder/profile.sscore.vars

# Combine log files chromosomes
touch $out_folder/plink.sscore.log
for chr in $(seq 1 22); do
  cat $out_folder/profile_chr$chr\.log >> $out_folder/plink.sscore.log
  cat $out_folder/flip_profile_chr$chr\.log >> $out_folder/plink.sscore.log
done
gzip $out_folder/plink.sscore.log

# remove all the chromosome specific files
rm $out_folder/*chr*

# gzip the final profile score as well
gzip $out_folder/profile.sscore

# remove the grs weights files
rm $out_folder/grs_weights_strand_flipped.txt.gz
rm $out_folder/grs_weights.txt.gz
