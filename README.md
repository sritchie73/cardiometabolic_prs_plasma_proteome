# README

This repository houses and documents the code used to generate the results in the study Ritchie SC *et al.* Integrative analysis of the plasma proteome and polygenic risk of cardiometabolic diseases. *bioRxiv*, doi: 10.1101/2019.12.14.876474 (https://www.biorxiv.org/content/10.1101/2019.12.14.876474v2.article-info).

This code has not been designed to regenerate the results as-is for third-parties. It has been written to run on a high-performance computing cluster at the University of Cambridge - it includes job submission scripts that are written specifically for this cluster's setup, and cannot which cannot be generalised. The underlying data is also not provided as part of this repo - these must be downloaded separately (see Underlying Data section below) and stored in the locations given in the hard-coded filepaths in these scripts. 

# Repository Structure

All scripts are housed under the src/ directory in this repository. At the top level, scripts are conceptually organised by analysis task, with one job submission script per task. For example, src/01_calculate_all_grs.sh is a job submission script that submits a sequence of jobs to the cluster, with the individual job scripts located in src/01_job_scripts/. These are designed to be part of a wider pipeline to test associations between an arbitrary number of polygenic scores (not included in this paper) with molecular measurements in INTERVAL from various high-throughput platforms (not included in this paper). Components of this pipeline that were not used in this paper are not included in this repositor, hence while the job scripts in this pipeline are sequential in run order, some numeric steps are not included in this repo. Some of these job dispatch scripts accept a "view file" as an argument (e.g. src/05_sensitivity_associations.sh), which contains a list of polygenic scores to analyse. The relevant view file needed for this paper is provided in the views/ folder, listing the five polygenic risk scores analysed.

# Software and versions used:

The following software and versions were used to run these scripts:

 - Scientific Linux release 7.7 (Nitrogen) (HPC operating system)
 - slurm version 19.05.5 (HPC queue manager and job submission system)
 - GNU bash version 4.2.46(2) (shell environment used to run bash scripts)
 - PLINK v1.90b6.10 64-bit (17 Jun 2019) (www.cog-genomics.org/plink/1.9/), aliased as plink1.9 in the scripts.
 - PLINK v2.00a2LM AVX2 Intel (24 Jul 2019) (www.cog-genomics.org/plink/2.0/), aliased as plink2 in the scripts.
 - R version 3.6, along with R packages:
     - data.table version 1.12.8
     - foreach version 1.4.4
     - doMC version 1.3.5
     - XML version 3.98-1.20
     - biomaRt version 2.40.3 (Bioconductor package)
     - openxlsx version 4.1.0.1
     - ggplot2 version 3.2.0
     - MendelianRandomization version 0.4.1
     - ggrepel version 0.8.1
 - The BGEN software suite (https://www.well.ox.ac.uk/~gav/bgen_format/software.html) including:
     - bgenix version 1.1.4
     - qctool version 2.0.5, alised as qctool2 in the scripts
     - ldstore version 1.1
 - SQLite version 3.30.1, aliased as sqlite3 in the scripts.

Inkscape version 0.92.3 was used to layout and annotate figures from the figure components generated within the R scripts. Microsoft Office Professional Plus 2016 was used to draft the manuscript (Microsoft Word) and curate supplemental tables (Microsoft Excel) on Windows 10 Enterprise edition.

# Underlying Data

All data used in this study is publicly available or deposited in a public repository. Genetic data, proteomic data, and basic cohort characteristics for the INTERVAL cohort are available via the European Genotype-phenome Archive (EGA) with study accession EGAS00001002555 (https://www.ebi.ac.uk/ega/studies/EGAS00001002555). Dataset access is subject to approval by a Data Access Committee: these data are not publicly available as they contain potentially identifying and sensitive patient information. All other data used in this study is publicly available without restriction. The PRS used for CAD in this study is available to download at the Polygenic Score Catalog with score accession PGS000018 (https://www.pgscatalog.org/pgs/PGS000018). The other PRSs used in this study are available through Figshare at https://dx.doi.org/10.6084/m9.figshare.11369103. GWAS summary statistics used to generate these PRSs are available to download through the GWAS Catalog (https://www.ebi.ac.uk/gwas/) with study accessions GCST008065 (CKD), GCST007517 (T2D), GCST006414 (Atrial fibrillation), and GCST006906 (Stroke, all causes). pQTL summary statistics used in Figure 3, Figure S4, Table S4, and Table S6 are available to download from https://www.phpc.cam.ac.uk/ceu/proteins/. GWAS summary statistics used for Mendelian randomisation (Figure 3, Table S6) are available to download through the GWAS Catalog (https://www.ebi.ac.uk/gwas/) with study accessions GCST004787 (CAD), GCST008065 (CKD), and GCST007518 (T2D). Correlations between tissue-specific expression and cardiometabolic traits in the hybrid mouse diversity panel used in Figure 4 and Table S10 are available to download through the systems genetics resource (https://systems.genetics.ucla.edu). Tissue-specific gene expression data and cardiometabolic phenotypes in the F2 cross of the inbred ApoE-/- C57BL/6J and C3H/HeJ strains used in Figure 4 and Table S10 are available to download from Sage BioNetworks at https://www.synapse.org/#!Synapse:syn4497. Tissue-specific gene expression data in C57BL/6J mice subject to dietary intervention generated for this study and used in Figure S5 and Table S11 is available in Table S12.
