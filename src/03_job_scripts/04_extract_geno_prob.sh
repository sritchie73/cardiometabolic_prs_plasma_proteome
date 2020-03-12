#!/bin/bash
chr=$SLURM_ARRAY_TASK_ID
omic_dir=analyses/omic_dirs/somalogic_proteins/

# Need to modify dispatch script if we want to run for other platforms
# due to runtime when loading bgen files

if [ -f $omic_dir/qtl_list.txt ] && [ ! -f $omic_dir/qtl_geno_prob.tsv ]; then
		# Create a recoded bim file
		Rscript src/03_job_scripts/04_helpers/bim_assign_blank_rsIDs.R \
			analyses/processed_genotypes/impute_chr$chr\_interval_filtered.bim \
			$omic_dir/chr$chr\.bim
		# symlink to rest of genotype data
		ln -s $(pwd)/analyses/processed_genotypes/impute_chr$chr\_interval_filtered.bed $omic_dir/chr$chr\.bed
		ln -s $(pwd)/analyses/processed_genotypes/impute_chr$chr\_interval_filtered.fam $omic_dir/chr$chr\.fam
		# extract variants
		plink2 --bgen data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/impute_$chr\_interval.bgen \
					 --sample data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/interval.samples \
           --memory 25000 --threads 1 \
					 --extract $omic_dir/qtl_variants.txt \
					 --export A --out $omic_dir/qtl_vars_chr$chr
		# extract alleles for each person
		plink1.9 --bfile $omic_dir/chr$chr \
						 --extract $omic_dir/qtl_variants.txt \
						 --recode --out $omic_dir/qtl_alleles_chr$chr
fi

