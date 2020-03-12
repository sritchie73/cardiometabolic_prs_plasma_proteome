#!/bin/bash

GRS=$1
Aptamer=$2
grs_folder=data/GRS_resources/$GRS
out_dir=analyses/grs_pqtl_removed/$GRS/$Aptamer

# Create a strand-flipped version of the grs
zcat $grs_folder/grs_weights.txt.gz |
  sed 's/A/V/g' | sed 's/T/X/g' | sed 's/C/Y/g' | sed 's/G/Z/g' |
  sed 's/V/T/g' | sed 's/X/A/g' | sed 's/Y/G/g' | sed 's/Z/C/g' \
  > $out_dir/grs_weights_strand_flipped.txt
gzip $out_dir/grs_weights_strand_flipped.txt

# Copy the GRS and edit both the GRS and strand flipped GRS to
# use common identifier format
cp $grs_folder/grs_weights.txt.gz $out_dir/
Rscript src/08_job_scripts/05_helpers/recode_grs_ids.R $out_dir/grs_weights.txt.gz
Rscript src/08_job_scripts/05_helpers/recode_grs_ids.R $out_dir/grs_weights_strand_flipped.txt.gz
