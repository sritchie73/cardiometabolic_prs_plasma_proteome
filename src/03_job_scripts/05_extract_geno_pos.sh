#!/bin/bash
n_var=$SLURM_ARRAY_TASK_ID
var=$(grep ":" analyses/processed_traits/somalogic_proteins/qtl_variants.txt | sed -n "${n_var}p")
pos=$(echo $var | cut -f 2 -d ":")
chr=$(echo $var | cut -f 1 -d ":" | sed "s/chr//")

# Need to modify dispatch script if we want to run for other platforms
# due to runtime when loading bgen files
omic_dir=analyses/processed_traits/somalogic_proteins

if [ -f $omic_dir/qtl_list.txt ] && [ ! -f $omic_dir/qtl_geno_prob.tsv ]; then
		# extract variants
		plink2 --bgen data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/impute_$chr\_interval.bgen \
					 --sample data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/interval.samples \
           --memory 25000 --threads 1 \
           --chr $chr --from-bp $pos --to-bp $pos \
					 --export A --out $omic_dir/qtl_vars_chr$chr\_$pos

    # Give column appropriate name
    Rscript src/03_job_scripts/05_helpers/fix_raw.R $omic_dir/qtl_vars_chr$chr\_$pos\.raw $chr $pos
fi

