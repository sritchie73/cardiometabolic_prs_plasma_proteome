#!/bin/bash

grs_name=$1
out_folder=analyses/GRS_profiles/$grs_name
geno_folder=analyses/processed_genotypes
grs_folder=data/GRS_resources/$grs_name

# Combine all the score files across chromosomes
Rscript src/01_job_scripts/04_helpers/chr_profile_combiner.R $grs_name

# Copy and combine plink log files from genotype qc stage
for chr in $(seq 1 22); do 
  cat $geno_folder/impute_chr$chr\_interval.log >> $out_folder/plink.qc.log
  cat $geno_folder/impute_chr$chr\_interval_deduped.log >> $out_folder/plink.qc.log
  cat $geno_folder/impute_chr$chr\_interval_filtered.log >> $out_folder/plink.qc.log
done
gzip $out_folder/plink.qc.log

cp $geno_folder/genotype_filters.txt $out_folder/genotype_filters.txt

# Summarise variant numbers
passed_qc=$(wc -l $geno_folder/impute_chr*_interval_filtered.bim | tail -1 | sed "s/^ \+//" | cut -f 1 -d " " | xargs printf "%'d")
echo "$passed_qc variants passed filters and qc." >> $out_folder/genotype_filters.txt

complement_snps=$(cat $out_folder/complementary_grs_snps_chr*.txt | paste -s -d + | bc | xargs printf "%'d")
echo "$complement_snps A/T or G/C variants found in GRS, these have been excluded." > $out_folder/grs_complement_snps.txt

for chr in $(seq 1 22); do 
  cat $out_folder/flip_profile_chr$chr.sscore.vars >> $out_folder/strand_mismatch.sscore.vars
done
gzip $out_folder/strand_mismatch.sscore.vars

# Combine all variant lists across chromosomes
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

# Remove chromosome specific files, unless they're the sscore files
# in which case we will save them for later
for ff in $out_folder/*chr*; do
  echo $ff | grep -q ".sscore$"
  if [[ $? -eq 0 ]]; then
    gzip $ff;
  else
    rm $ff
  fi
done

# gzip the final profile score as well
gzip $out_folder/profile.sscore

# remove the grs weights files with recoded ids
rm $out_folder/grs_weights_strand_flipped.txt.gz
rm $out_folder/grs_weights.txt.gz
