#!/bin/bash

GRS=$1
Aptamer=$2
chr=$SLURM_ARRAY_TASK_ID

out_dir=analyses/grs_pqtl_removed/$GRS/$Aptamer/

# Check if this chromosome is one that has variants to be removed
pqtl_chr=$(cat $out_dir/chrs_to_filter.txt)
chr_has_pqtl=false
for elem in $pqtl_chr; do
  if [ $elem -eq $chr ]; then
    chr_has_pqtl=true
  fi
done

# Filter out pQTLs and tagging variants
if $chr_has_pqtl; then
	plink1.9 --bfile analyses/grs_pqtl_removed/chr$chr \
           --exclude $out_dir/pqtl_tags.txt \
					 --make-bed --out $out_dir/chr$chr

	# Cleanup
	rm $out_dir/chr$chr\.nosex
	rm $out_dir/chr$chr\.log

	# Revert bim file variant coding so that GRS variants will be match
	Rscript src/08_job_scripts/04_helpers/recode_bim.R $out_dir/chr$chr\.bim
fi

