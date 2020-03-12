#!/bin/bash

if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then
  echo "Script must be submitted using sbatch with the --array argument"
  exit 1
fi
chr=$SLURM_ARRAY_TASK_ID

# file paths
out_folder="analyses/processed_genotypes"
snpstats_folder="data/INTERVAL/post_qc_data/reference_files/genetic/reference_files/"

# Customisable parameters. Comment out options you don't want to use
memory_limit=$1
dosage_certainty=0.1 # plink makes hard calls based on genotype certainty, or sets to missing based on this threshold 
info=0.3
# maf=0.01
# geno=0.1 # genotype missingness filter

# The first step is to convert the BGEN files to .bed, .bim, and .bam files.
# This means future plink calls do not need to do the conversion internally,
# and allows us to subsequently process or recode variant IDs trivially.
# Plink calls are built using arrays where the script needs
# to dynamically detect arguments where plink would error on blank input.
convert_cmd[0]='plink2'
convert_cmd[1]='--bgen data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/impute_$chr\_interval.bgen'
convert_cmd[2]='--sample data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/interval.samples'
convert_cmd[3]='--memory $memory_limit --threads 1'
if [ ! -z ${dosage_certainty+x} ]; then # if not set, plink will use 0.1 anyway
  convert_cmd[4]='--import-dosage-certainty $dosage_certainty'
fi
convert_cmd[5]='--make-bed --out $out_folder/impute_chr$chr\_interval'

eval ${convert_cmd[@]}

# Give each variant a strictly unique id so that we can later identify and remove true
# duplicates (i.e. variants sharing a chromosome, position, and allele combinations)
Rscript src/01_job_scripts/01_helpers/make_ids_unique.R $out_folder/impute_chr$chr\_interval.bim

# Create an additional file containing the INFO scores for each variant. We
# will use this to identify and resolve duplicates, as well as to filter by
# INFO score if needed. 
cut -f 2 -d " " $out_folder/impute_chr$chr\_interval.bim > $out_folder/ids_chr$chr.txt # uids
tail -n +2 $snpstats_folder/impute_$chr\_interval.snpstats | cut -f 3-6,19 > $out_folder/info_chr$chr.txt
paste $out_folder/ids_chr$chr.txt $out_folder/info_chr$chr.txt > $out_folder/info_filter_chr$chr.txt
rm $out_folder/ids_chr$chr.txt $out_folder/info_chr$chr.txt

# Resolve any duplicate variants by keeping only the one with the highest INFO 
# score. Variants are matched based on chromosome, position, and alleles for the
# purposes of identifying duplicates, rather than by rsID, because we may want to
# convert to this format later and don't want to lose variants without rsIDs in
# INTERVAL (currently given an rsID of ".").
Rscript src/01_job_scripts/01_helpers/flag_duplicates_ids_to_remove.R \
        $out_folder/info_filter_chr$chr.txt \
        $out_folder/duplicates_to_remove_chr$chr.txt

plink2 --bfile $out_folder/impute_chr$chr\_interval \
         --memory $memory_limit --threads 1 \
         --exclude $out_folder/duplicates_to_remove_chr$chr.txt \
         --make-bed --out $out_folder/impute_chr$chr\_interval_deduped

# remove previous larger initial bed/bim/fam files to save disk space
rm $out_folder/impute_chr$chr\_interval.bed
rm $out_folder/impute_chr$chr\_interval.bim
rm $out_folder/impute_chr$chr\_interval.fam

# Now apply the filtering. Built with command array to allow for varying of
# which filters get applied
filter_cmd[0]='plink1.9' # plink1.9 required for --qual-scores
filter_cmd[1]='--bfile $out_folder/impute_chr$chr\_interval_deduped'
filter_cmd[2]='--memory $memory_limit --threads 1'
if [ ! -z ${info+x} ]; then
  filter_cmd[3]='--qual-scores $out_folder/info_filter_chr$chr.txt 6 1'
  filter_cmd[4]='--qual-threshold $info' # INFO score threshold
fi
if [ ! -z ${geno+x} ]; then
  filter_cmd[5]='--geno $geno' # genotype missingness
fi
if [ ! -z ${maf+x} ]; then
  filter_cmd[6]='--maf $maf'
fi
filter_cmd[7]='--make-bed --out $out_folder/impute_chr$chr\_interval_filtered'

eval ${filter_cmd[@]}

# remove previous larger deduped bed/bim/fam files to save disk space
rm $out_folder/impute_chr$chr\_interval_deduped.bed
rm $out_folder/impute_chr$chr\_interval_deduped.bim
rm $out_folder/impute_chr$chr\_interval_deduped.fam

# Remove suffix used to identify and remove duplicates from the variant IDs in the .bim file:
sed -i "s/:[0-9]*\t/\t/" $out_folder/impute_chr$chr\_interval_filtered.bim

# Output file containing filters in effect:
if [ $chr -eq 22 ]; then 
  echo "Genotype filtering in effect:" >> $out_folder/genotype_filters.txt
  
  if [ -z ${dosage_certainty+x} ]; then
    dosage_certainty=0.1
  fi 
  echo "Hard call made on genotype based on allele with highest certainty. Set to missing where this is < $dosage_certainty." >>  $out_folder/genotype_filters.txt
  echo "Duplicate variants resolved by keeping only the one with the highest INFO." >> $out_folder/genotype_filters.txt

  if [ ! -z ${geno+x} ]; then
    echo "Variants excluded where genotype missingness < $geno."
  fi
  
  if [ ! -z ${info+x} ]; then
    echo "Variants excluded where INFO < $info." >>  $out_folder/genotype_filters.txt
  fi

  if [ ! -z ${maf+x} ]; then
    echo "Variants excluded where their MAF < $maf." >> $out_folder/genotype_filters.txt
  fi
fi
