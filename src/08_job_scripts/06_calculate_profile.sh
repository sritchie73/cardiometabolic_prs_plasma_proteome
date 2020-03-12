#!/bin/bash

if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then
  echo "Script must be submitted using sbatch with the --array argument"
  exit 1
fi

grs_name=$1
Aptamer=$2
memory_limit=$3
chr=$SLURM_ARRAY_TASK_ID

profile_folder="analyses/GRS_profiles/$grs_name"
out_folder="analyses/grs_pqtl_removed/$grs_name/$Aptamer"

# First, check if this is one of the chromosomes where variants have been removed.
# If not, we can just symlink to the original chromosome score file to save time!
pqtl_chr=$(cat $out_folder/chrs_to_filter.txt)
done_already=true
for chrIdx in $pqtl_chr; do
  if [ $chrIdx -eq $chr ]; then
    done_already=false
  fi
done
if $done_already; then
  # Symlink to score file(s) for this chromosome
  profile_file=$profile_folder/profile_chr$chr\.sscore.gz
  flip_file=$profile_folder/flip_profile_chr$chr\.sscore.gz
	ln -s $PWD/$profile_file $out_folder/$(basename $profile_file)
  if [ -f $flip_file ]; then
    ln -s $PWD/$flip_file $out_folder/$(basename $flip_file)
    has_chr_score=true
  fi

  # Also get the sscore components so we can count
  # the number of variants before and after pQTL removal
  zcat $profile_folder/profile.sscore.vars.gz | grep "^${chr}:" > $out_folder/profile_chr$chr\.sscore.vars
  exit 0
fi

# First we need to convert all identifiers to chromosome:position:allele1:allele2
# for matching since there can be differences in how dataset code ids for variants
# with no rsID, or where the alleles themselves are part of the identifiers. Here,
# alleles are alphabetically sorted for guaranteed matching.

# First edit the BIM files directly to code the identifiers. We can then symlink the
# to the .bed and .fam files rather than copying them
Rscript src/08_job_scripts/06_helpers/replace_variant_ids_bim.R "$out_folder/chr${chr}.bim"

# Identify A/T and G/C variants, which will be ignored when calculating the GRS
Rscript src/08_job_scripts/06_helpers/identify_complementary_snps.R \
        "$out_folder/chr${chr}.bim" \
        "$out_folder/complementary_snps_chr${chr}.txt"

# Identify how many of the complementary alleles in INTERVAL are in the GRS
# and log for later collation
Rscript src/08_job_scripts/06_helpers/identify_GRS_complementary_snps.R \
        "$out_folder/complementary_snps_chr${chr}.txt" \
        "$out_folder/grs_weights.txt.gz" \
        "$out_folder/complementary_grs_snps_chr${chr}.txt"

# Build command array for calculating the GRS score levels. We will do this twice:
# the second time with alleles in the GRS strand flipped to capture any variants
# that were on the opposite strand in the dataset(s) used to generate the GRS
cmd[0]='plink2'
cmd[1]='--bfile $out_folder/chr$chr'
cmd[2]='--memory $memory_limit --threads $SLURM_CPUS_ON_NODE'
cmd[3]='--exclude $out_folder/complementary_snps_chr$chr.txt' # works ok if file is present but empty
cmd[4]='--score $weights_file 1 4 6 ignore-dup-ids list-variants header cols=nmissallele,scoresums'
cmd[5]='--out $out_folder/$out_file'

# Run once using the weights as is:
weights_file="$out_folder/grs_weights.txt.gz"
out_file="profile_chr$chr"
eval ${cmd[@]}

# Same again, just for opposite strand in cases there are strand flip
# issues for any variants. We must exclude any GRS variants successfully
# matched in the original strand orientation to prevent cancelling out
# alleles were the reference and alternative are complementary bases
cat $out_folder/$out_file.sscore.vars >> $out_folder/complementary_snps_chr$chr.txt
weights_file="$out_folder/grs_weights_strand_flipped.txt.gz"
out_file="flip_profile_chr$chr"

eval ${cmd[@]}

# Remove links to bed/bim/fam files to clear up directory
rm "$out_folder/chr${chr}.bed"
rm "$out_folder/chr${chr}.bim"
rm "$out_folder/chr${chr}.fam"

