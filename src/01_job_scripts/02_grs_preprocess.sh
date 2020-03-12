#!/bin/bash

grs_name=$1
grs_folder=data/GRS_resources/$grs_name
out_folder=analyses/GRS_profiles/$grs_name

# Create a strand-flipped version of the grs
zcat $grs_folder/grs_weights.txt |
  sed 's/A/V/g' | sed 's/T/X/g' | sed 's/C/Y/g' | sed 's/G/Z/g' |
  sed 's/V/T/g' | sed 's/X/A/g' | sed 's/Y/G/g' | sed 's/Z/C/g' \
  > $out_folder/grs_weights_strand_flipped.txt
gzip $out_folder/grs_weights_strand_flipped.txt

# Copy the GRS and edit both the GRS and strand flipped GRS to
# use common identifier format
cp $grs_folder/grs_weights.txt.gz $out_folder/
Rscript src/01_job_scripts/02_helpers/recode_grs_ids.R $out_folder/grs_weights.txt.gz
Rscript src/01_job_scripts/02_helpers/recode_grs_ids.R $out_folder/grs_weights_strand_flipped.txt.gz

