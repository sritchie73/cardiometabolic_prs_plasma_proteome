library(data.table)
library(foreach)
library(openxlsx)
source("src/utilities/flip_strand.R")
source("src/utilities/prot_pval.R")

# Determine GRS we're working with
if (!exists("grs")) {
  view_file <- commandArgs(trailingOnly=TRUE)[1]
  view <- fread(view_file)
  task <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  grs <- view[task, GRS_name]
}

# Create directory to hold collated results
out_dir <- sprintf("analyses/mendelian_randomisation/%s", grs)

# Resource directories:
pqtl_dir <- "data//full_pQTL_summary_stats"
gwas_dir <- sprintf("analyses/mendelian_randomisation/%s/gwas_summary_stats", grs)
GWASs <- list.files(path=gwas_dir)
snp_dir <- "data/INTERVAL/reference_files/imputed_genotypes"
pqtl_ld <- "analyses/mendelian_randomisation/pqtls/conditional_pQTLs_tag_LD"
cis_ld <- "analyses/mendelian_randomisation/pqtls/SOMAMERs_cis_LD"

# Determine aptamers we're working with
info <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
assocs <- fread(sprintf("analyses/univariate_associations/%s/somalogic_proteins/associations.tsv", grs))
assocs <- assocs[info, on=.(trait=variable), nomatch=0]
assocs <- assocs[Type == "Protein"]
prot_assocs <- assocs[, .(pval = prot_pvalue(pval, beta)),by=.(Target, Gene.Name, UniProt=UniProt.Id.Current.at.Uniprot)]
prot_assocs[, fdr := p.adjust(pval, method="fdr")]
#prot_assocs <- prot_assocs[fdr < 0.05]
assocs <- assocs[prot_assocs[,.(Target, Gene.Name, UniProt)], 
                 on = .(Target, Gene.Name, UniProt.Id.Current.at.Uniprot=UniProt)]

# Check if there is anything to run:
if (nrow(assocs) == 0) {
  quit(save="no", status=0)
}
if (file.exists(sprintf("%s/instruments.tsv", out_dir))) {
  quit(save="no", status=0)
}

# Load published pQTL information
pQTLs <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=5, startRow=3)
pQTLs <- as.data.table(pQTLs)
pQTLs <- pQTLs[SOMAmer.ID %in% assocs$SOMAMER_ID]

# Add information about whether each variant is cis or trans
head_row <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, rows=5:6, fillMergedCells=TRUE)
head_row <- gsub(".NA$", "", paste(colnames(head_row), as.vector(head_row[1,]), sep="."))
pQTL_info <- read.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=4, startRow=6)
colnames(pQTL_info) <- head_row
pQTL_info <- as.data.table(pQTL_info)

pQTLs[pQTL_info, on = .(SOMAmer.ID, Sentinel.variant=`Sentinel.variant*`), type := `cis/.trans`]

# Load sub-threshold cis-pQTL information
cis_sub <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
cis_sub <- cis_sub[SOMAMER_ID %in% assocs$SOMAMER_ID]

# Setup - iterate per chromosome for efficiency:
rbindf <- function(...) { rbindlist(list(...), fill=TRUE, use.names=TRUE) }
last_chr <- -1
snpstats <- data.table()

# Determine chromosomes
chrs <- sort(unique(c(cis_sub$chr, pQTLs$Conditional.variant.Chr)))

# Check whether any proteins have cis-pQTLs
if (length(chrs) == 0) {
  quit(save="no", status=0)
}

# Load snpstats for INTERVAL
chr_snpstats <- foreach(chrIdx = chrs, .combine=rbind) %do% {
  snpstats <- fread(sprintf("%s/impute_%s_interval.snpstats", snp_dir, chrIdx))
  snpstats[, chromosome := as.character(chromosome)]
  snpstats
}

# Map pQTL summary statistics to GWAS summary statistics
instruments <- foreach(gwas = GWASs, .combine=rbindf) %do% {
  # Load GWAS summary stats
  gwas_chr_stats <- foreach(chrIdx = chrs, .combine=rbind) %do% {
    gwas_stats <- fread(sprintf("%s/%s/chr%s.txt.gz", gwas_dir, gwas, chrIdx))
    gwas_stats <- gwas_stats[, .(chr, pos, EA, OA, EAF=MAF, effect, se, P)]
    gwas_stats[, chr := as.character(chr)]
    gwas_stats
  } 
  
  # Filter out variants that can't be reliably mapped due to strand ambiguity
  gwas_chr_stats <-  gwas_chr_stats[EA != flip_strand(OA)]

  # Map IVs for each protein
  prots <- unique(assocs[,.(Target, UniProt.Id.Current.at.Uniprot, Gene.Name)])
  foreach(protIdx = prots[,.I], .combine=rbindf) %do% {
    prot <- prots[protIdx]
    # First, do published pQTLs for any aptamers for that protein.
    prot_aptamers <- assocs[prot, on = .(Target, UniProt.Id.Current.at.Uniprot, Gene.Name), unique(SOMAMER_ID)]
    published_pQTLs <- pQTLs[SOMAmer.ID %in% prot_aptamers]
    pqtl_locs <- unique(published_pQTLs[,.(chr=Conditional.variant.Chr, pos=Conditional.variant.Pos)])
    pub_map <- foreach(pqtl_loc_idx = pqtl_locs[,.I], .combine=rbindf) %do% {
      pqtl_loc <- pqtl_locs[pqtl_loc_idx]
      # Collate variants in LD
      ld <- fread(sprintf("%s/chr%s/%s.ld", pqtl_ld, pqtl_loc$chr, pqtl_loc$pos)) 
      ld <- ld[, .(pos = position2, r2)]
      ld <- rbind(ld, data.table(pos=pqtl_loc$pos, r2=1))

      # Load pQTL summary stats for each aptamer for the chromosome this pQTL is on
      pqtl_chr_stats <- foreach(apt = prot_aptamers, .combine=rbindf) %do% {
         apt_pqtl_stats <- fread(sprintf("%s/%s/%s_chrom_%s_meta_1.tbl.gz", pqtl_dir, apt, apt, pqtl_loc$chr))
         apt_pqtl_stats[, SOMAMER_ID := apt]
         apt_pqtl_stats
      }
      pqtl_chr_stats <- pqtl_chr_stats[, .(SOMAMER_ID, chr=as.character(chromosome), pos=position, 
                                           EA=toupper(Allele1), OA=toupper(Allele2),
                                           effect=Effect, se=StdErr, P=10^`log(P)`)]

      # Filter to variants in LD with the pQTL 
      pqtl_chr_stats <- pqtl_chr_stats[ld, on=.(pos)]

      # Add effect size allele frequencies and RSIDs
      pqtl_chr_stats[chr_snpstats, on = .(chr=chromosome, pos=position, EA=minor_allele, OA=major_allele), c("EAF", "RSID") := .(MAF, i.RSID)]
      pqtl_chr_stats[chr_snpstats, on = .(chr=chromosome, pos=position, EA=major_allele, OA=minor_allele), c("EAF", "RSID") := .(1 - MAF, i.RSID)]

      # Get pQTLs and proxies that can be directly matched between datasets:
      strand_ok_col_ok <- merge(pqtl_chr_stats, gwas_chr_stats, 
                                by=c("chr", "pos", "EA", "OA"), suffixes=c(".pqtl", ".gwas"))

      # Handle cases where the effect allele is the other allele in the GWAS dataset
      strand_ok_col_swap <- merge(pqtl_chr_stats, gwas_chr_stats, 
                                  by.x=c("chr", "pos", "EA", "OA"), by.y=c("chr", "pos", "OA", "EA"), 
                                  suffixes=c(".pqtl", ".gwas"))
      strand_ok_col_swap[, c("effect.gwas", "EAF.gwas") := .(-effect.gwas, 1 - EAF.gwas)]

      # Now get ones that can be directly matched after flipping the strand:
      pqtl_chr_stats[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]
      strand_flip_col_ok <- merge(pqtl_chr_stats, gwas_chr_stats, 
                                  by=c("chr", "pos", "EA", "OA"), suffixes=c(".pqtl", ".gwas"))
      strand_flip_col_ok[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

      # Strand flip, effect allele = other allele
      strand_flip_col_swap <- merge(pqtl_chr_stats, gwas_chr_stats, 
                                    by.x=c("chr", "pos", "EA", "OA"), by.y=c("chr", "pos", "OA", "EA"), 
                                    suffixes=c(".pqtl", ".gwas"))
      strand_flip_col_swap[, c("effect.gwas", "EAF.gwas") := .(-effect.gwas, 1 - EAF.gwas)]
      strand_flip_col_swap[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

      # Reset pqtl_chr_stats to original strand
      pqtl_chr_stats[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

      # Alleles always oriented to effect allele and strand of pQTL study:
      combined <- rbind(strand_ok_col_ok, strand_ok_col_swap, 
                        strand_flip_col_ok, strand_flip_col_swap)

      # Add pQTL summary stats for variants with no-match in the GWAS
      pqtl_only <- pqtl_chr_stats[!combined, on=.(pos)]
      setnames(pqtl_only, 
               c("effect", "se", "P", "EAF"), 
               c("effect.pqtl", "se.pqtl", "P.pqtl", "EAF.pqtl"))
      combined <- rbind(combined, pqtl_only, fill=TRUE)

      # Collate information about pQTL and matches:
      if (combined[pos != pqtl_loc$pos, .N] > 0) {
        match_stats <- combined[pos != pqtl_loc$pos, 
           .(proxies = .N, matched = sum(!is.na(P.gwas))),
           by=SOMAMER_ID] 
      } else {
        match_stats <- data.table(SOMAMER_ID=prot_aptamers, proxies=0, matched=NA_integer_)
      }

      # Average variant effects across aptamers - note this assumes the variant has the
      # same effect allele across (all) aptamers - the end user should double check this!
      combined_prot <- combined[, .(effect.pqtl=mean(effect.pqtl), se.pqtl=mean(se.pqtl),
                                   P.pqtl=prot_pvalue(P.pqtl, effect.pqtl)),
                                   by=.(chr, pos, EA, OA, r2, EAF.pqtl, RSID, EAF.gwas,
                                   effect.gwas, se.gwas, P.gwas)]

      # Take the best variant (pQTL, or if not matched highest LD)
      if (nrow(combined_prot[!is.na(effect.gwas)]) > 0) {
        combined_prot <- combined_prot[!is.na(effect.gwas)]
        combined_prot <- combined_prot[order(abs(pqtl_loc$pos-pos))][order(-r2)]
        combined_prot <- combined_prot[1]
      } else {
        combined_prot <- combined_prot[pos == pqtl_loc$pos]
      }

      # Merge in aptamer level information
      setnames(combined_prot, c("effect.pqtl", "se.pqtl", "P.pqtl"),  
               c("Prot.Effect.pQTL", "Prot.SE.pQTL", "Prot.P.pQTL"))
      setnames(combined, c("effect.pqtl", "se.pqtl", "P.pqtl"),
               c("Apt.Effect.pQTL", "Apt.SE.pQTL", "Apt.P.pQTL"))
      combined <- combined_prot[combined, on = .(RSID, chr, pos, EA, OA, EAF.pqtl, r2, 
                                                 EAF.gwas, effect.gwas, se.gwas, P.gwas),
                                nomatch=0, allow.cartesian=TRUE]

      # Add match stats
      combined <- cbind(combined, data.table("match"=NA_character_))
      combined <- combined[match_stats, on = .(SOMAMER_ID)]

      # Which aptamer is this a pQTL for?
      combined[, aptamer_pqtl := FALSE]
      pqtl_loc[, chr := as.numeric(chr)]
      pqtl_apts <- published_pQTLs[pqtl_loc, on = .(Conditional.variant.Chr=chr, Conditional.variant.Pos=pos), SOMAmer.ID]
      combined[SOMAMER_ID %in% pqtl_apts, aptamer_pqtl := TRUE]
      pqtl_loc[, chr := as.character(chr)]

      # Determine match or non-match conditions
      combined[pqtl_loc$pos == pos & !is.na(P.gwas), match := "direct"]
      combined[pqtl_loc$pos != pos & !is.na(P.gwas), match := "proxy"]
      combined[pqtl_loc$pos == pos & is.na(P.gwas), match := "no matches"]
      combined[match == "direct", r2 := NA_real_]

      # Add additional annotation information
      anno <- published_pQTLs[Conditional.variant.Pos == pqtl_loc$pos, 
                              .(proxy_for=Conditional.variant, sentinel=Sentinel.variant,
                                type, type2="published")]
      anno[proxy_for == sentinel, sentinel := NA_character_]
      combined <- cbind(combined, anno)
      combined[pos == pqtl_loc$pos, proxy_for := NA_character_]
      combined[, chr := pqtl_loc$chr]
      combined[, c("Target", "UniProt", "Gene") := prot]
      combined[, GWAS := gwas]
      combined[, PRS := grs]
      return(combined)
    }

    if (is.null(pub_map)) {
      pub_map <- data.table(PRS = character(0), GWAS = character(0), Target = character(0),
                            UniProt = character(0), Gene = character(0),
    											  RSID = character(0), chr = integer(0), pos = integer(0), 
                            EA = character(0), OA = character(0), type = character(0),
                            type2 = character(0), sentinel = character(0), match = character(0),
                            proxy_for = character(0), r2 = numeric(0), proxies = integer(0),
                            matched = integer(0), EAF.pqtl = numeric(0),
                            Prot.Effect.pQTL = numeric(0), Prot.SE.pQTL = numeric(0), Prot.P.pQTL = numeric(0),
                            Apt.Effect.pQTL = numeric(0), Apt.SE.pQTL = numeric(0), Apt.P.pQTL = numeric(0),
                            EAF.gwas = logical(0), effect.gwas = numeric(0), se.gwas = numeric(0),
                            Prot.Effect.pQTL = numeric(0), Prot.SE.pQTL = numeric(0), Prot.P.pQTL = numeric(0),
                            P.gwas = numeric(0), SOMAMER_ID = character(0), aptamer_pqtl = logical(0))
    }

    # Now handle cis-pQTLs for the aptamers for this protein
    prot_aptamers <- assocs[prot, on = .(Target, UniProt.Id.Current.at.Uniprot, Gene.Name), unique(SOMAMER_ID)]
    prot[info, on = .(Target, UniProt.Id.Current.at.Uniprot, Gene.Name), c("chr", "start") := .(i.chr, i.start)]
    prot <- prot[,.(Gene.Name = strsplit(Gene.Name, "\\|")[[1]], 
                    chr = strsplit(chr, "\\|")[[1]],
                    start = strsplit(start, "\\|")[[1]]),
                 by = .(Target, UniProt.Id.Current.at.Uniprot)]
    prot <- prot[, .(Gene.Name = strsplit(Gene.Name, "\\;")[[1]],
                     chr = strsplit(chr, "\\;")[[1]],
                     start = strsplit(start, "\\|")[[1]]),
                 by = .(Target, UniProt.Id.Current.at.Uniprot)]
    cis_map <- foreach(gene_idx = prot[,.I], .combine=rbindf) %do% {
      gene <- prot[gene_idx]
      # Obtain cis pQTL summary statistics for all aptamers, filtering to cis-pQTLs significant
      # for both aptamers
      pqtl_chr_stats <- cis_sub[SOMAMER_ID %in% prot_aptamers]
      all_apt <- pqtl_chr_stats[, .N, by=.(chr, pos, EA, OA)][N == length(prot_aptamers)]
      pqtl_chr_stats <- pqtl_chr_stats[all_apt, on = .(chr, pos, EA, OA)]
      pqtl_chr_stats[, chr := as.character(chr)]

			# Go to next protein if no cis-PQTLs
			if (nrow(pqtl_chr_stats) == 0) return(NULL)

			# Add effect size allele frequencies and RSIDs
			pqtl_chr_stats[chr_snpstats, on = .(chr=chromosome, pos=position, EA=minor_allele, OA=major_allele), c("EAF", "RSID") := .(MAF, i.RSID)]
			pqtl_chr_stats[chr_snpstats, on = .(chr=chromosome, pos=position, EA=major_allele, OA=minor_allele), c("EAF", "RSID") := .(1 - MAF, i.RSID)]

			# Filter out variants with complementary alleles - we didnt 
      # calculate LD for these
      pqtl_chr_stats <- pqtl_chr_stats[EA != flip_strand(OA)]

			# Get pQTLs and proxies that can be directly matched between datasets:
			strand_ok_col_ok <- merge(pqtl_chr_stats, gwas_chr_stats, 
																by=c("chr", "pos", "EA", "OA"), suffixes=c(".pqtl", ".gwas"))

			# Handle cases where the effect allele is the other allele in the GWAS dataset
			strand_ok_col_swap <- merge(pqtl_chr_stats, gwas_chr_stats, 
																	by.x=c("chr", "pos", "EA", "OA"), by.y=c("chr", "pos", "OA", "EA"), 
																	suffixes=c(".pqtl", ".gwas"))
			strand_ok_col_swap[, c("effect.gwas", "EAF.gwas") := .(-effect.gwas, 1 - EAF.gwas)]
		 
			# Now get ones that can be directly matched after flipping the strand:
			pqtl_chr_stats[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]
			strand_flip_col_ok <- merge(pqtl_chr_stats, gwas_chr_stats, 
																	by=c("chr", "pos", "EA", "OA"), suffixes=c(".pqtl", ".gwas"))
			strand_flip_col_ok[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

			# Strand flip, effect allele = other allele
			strand_flip_col_swap <- merge(pqtl_chr_stats, gwas_chr_stats, 
																		by.x=c("chr", "pos", "EA", "OA"), by.y=c("chr", "pos", "OA", "EA"), 
																		suffixes=c(".pqtl", ".gwas"))
			strand_flip_col_swap[, c("effect.gwas", "EAF.gwas") := .(-effect.gwas, 1 - EAF.gwas)]
			strand_flip_col_swap[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

			# Reset pqtl_chr_stats to original strand
			pqtl_chr_stats[, c("OA", "EA") := .(flip_strand(OA), flip_strand(EA))]

			# Alleles always oriented to effect allele and strand of pQTL study:
			combined <- rbind(strand_ok_col_ok, strand_ok_col_swap, 
												strand_flip_col_ok, strand_flip_col_swap)

      # Go to next protein if no cis-pQTLs could be mapped
      if (nrow(combined) == 0) return(NULL)

      # create averages across aptamers
      combined_prot <- combined[, .(effect.pqtl=mean(effect.pqtl), se.pqtl=mean(se.pqtl),
                                     P.pqtl=prot_pvalue(P.pqtl, effect.pqtl)),
                                   by=.(chr, pos, EA, OA, EAF.pqtl, RSID, EAF.gwas,
                                        effect.gwas, se.gwas, P.gwas)]

			# Load in LD-map for this aptamer
      ld <- foreach(aptamer = prot_aptamers, .combine=rbind) %do% {
         fread(sprintf("%s/%s.ld", cis_ld, aptamer)) 
      }
      ld <- unique(ld)
			ld[, r2 := correlation^2]

			# Drop any variants in LD with a published pQTL
      if (!is.null(pub_map)) {
				cis_pub <- pub_map[chr == gene$chr & 
                           pos > as.numeric(gene$start) - 1e6 & 
                           pos < as.numeric(gene$start) + 1e6]
				if (nrow(cis_pub) > 0) {
					for (rIdx in 1:nrow(cis_pub)) {
						pqtl_pos <- cis_pub[rIdx, pos]
						pos_in_ld <- rbind(
							ld[position1 == pqtl_pos & r2 > 0.1, .(pos=position2)],
							ld[position2 == pqtl_pos & r2 > 0.1, .(pos=position1)],
							cis_pub[rIdx, .(pos)])
						pos_in_ld <- unique(pos_in_ld)
						combined_prot <- combined_prot[!pos_in_ld, on = .(pos)]
					}
				}
      }

      # Go to next aptamer if no cis-pQTLs left
      if (nrow(combined_prot) == 0) return(NULL)

			# Run a step forward procedure in which we take the top remaining pQTL,
			# then exclude variants in LD
			combined_prot <- combined_prot[order(-abs(effect.pqtl))][order(P.pqtl)]
			sub_threshold <- combined_prot[0]
			while (nrow(combined_prot) > 0) {
				top_pqtl <- combined_prot[1]
				sub_threshold <- rbind(sub_threshold, top_pqtl) 
				pos_in_ld <- rbind(
					ld[position1 == top_pqtl$pos & r2 > 0.1, .(pos=position2)],
					ld[position2 == top_pqtl$pos & r2 > 0.1, .(pos=position1)],
					top_pqtl[,.(pos)])
				pos_in_ld <- unique(pos_in_ld)
				combined_prot <- combined_prot[!pos_in_ld, on = .(pos)]
			}

      # Merge in aptamer level information
      setnames(sub_threshold, c("effect.pqtl", "se.pqtl", "P.pqtl"),  
               c("Prot.Effect.pQTL", "Prot.SE.pQTL", "Prot.P.pQTL"))
      setnames(combined, c("effect.pqtl", "se.pqtl", "P.pqtl"),
               c("Apt.Effect.pQTL", "Apt.SE.pQTL", "Apt.P.pQTL"))
      sub_threshold <- sub_threshold[combined, on = .(RSID, chr, pos, EA, OA, EAF.pqtl, 
                                                     EAF.gwas, effect.gwas, se.gwas, P.gwas),
                                               nomatch=0, allow.cartesian=TRUE]

			# Add annotation information
			sub_threshold[, type := "cis"]
			sub_threshold[, type2 := "sub-threshold"]
			sub_threshold[, GWAS := gwas]
			sub_threshold[, chr := gene$chr]
			sub_threshold[, PRS := grs]
      sub_threshold[, c("Target", "UniProt", "Gene") := gene[,.(Target, UniProt.Id.Current.at.Uniprot, Gene.Name)]]

      return(sub_threshold)
    }
    return(rbindf(pub_map, cis_map))
  }
}

# Orient effect sizes to minor allele in INTERVAL
instruments[EAF.pqtl > 0.5, 
  c("EA", "OA", "Prot.Effect.pQTL", "Apt.Effect.pQTL", "effect.gwas", "EAF.pqtl", "EAF.gwas") :=
  .(OA, EA, -Prot.Effect.pQTL, -Apt.Effect.pQTL, -effect.gwas, 1 - EAF.pqtl, 1 - EAF.gwas)]

# Add alternative identifiers
instruments[, ALT_ID := paste0("chr", chr, ":", pos, ":", EA, ":", OA)]
instruments[is.na(RSID), RSID := "."]

# Impose some sanity on the column and row order:
instruments <- instruments[, .(
  PRS, GWAS, Target, UniProt, Gene, RSID, ALT_ID, chr, pos, EA, OA, 
  type, type2, sentinel,
  match, proxy_for, r2_with_proxy=r2,
  proxies, matched,
  EAF.pqtl, Prot.Effect.pQTL, Prot.SE.pQTL, Prot.P.pQTL,
  EAF.gwas, effect.gwas, se.gwas, P.gwas,
  SOMAMER_ID, Apt.Effect.pQTL, Apt.SE.pQTL, Apt.P.pQTL, aptamer_pqtl)]
instruments <- instruments[order(type)][order(-type2)][order(Gene)][order(GWAS)][order(PRS)]

# Write out
fwrite(instruments, sep="\t", quote=FALSE, file=sprintf("%s/instruments.tsv", out_dir))

# Collate information about LD between instruments
out_dir <- sprintf("%s/LD", out_dir)
foreach(gwas = unique(instruments$GWAS), .errorhandling="remove") %do% {
  prots <- unique(instruments[,.(Target, UniProt, Gene)])
  foreach(protIdx = prots[,.I], .errorhandling="remove") %do% {
    prot <- prots[protIdx]
    sub_dt <- instruments[GWAS == gwas][prot, on = .(Target, UniProt, Gene)]

    # Create output directory
    out_name <- prot[,paste(Target, UniProt, Gene, sep="_")]
    out_name <- gsub(";", ".", out_name)
    out_name <- gsub("\\|", ".", out_name)
    out_name <- gsub(" ", ".", out_name)
    out_name <- gsub(",", ".", out_name)
    ld_out <- sprintf("%s/%s/%s", out_dir, out_name, gwas)
    dir.create(ld_out, recursive=TRUE, showWarnings=FALSE)

    # Filter to unique variants
    sub_dt <- unique(sub_dt, by=c("chr", "pos", "EA", "OA"))

    #preallocate ld-matrix:
    sub_dt <- sub_dt[order(pos)][order(as.numeric(chr))]
    ld_mat <- matrix(0, nrow=nrow(sub_dt), ncol=nrow(sub_dt), 
                     dimnames=list(sub_dt$ALT_ID, sub_dt$ALT_ID))
    
    # cis-LD is most straightforward, we already have that information calculated:
    cis_dt <- sub_dt[type == "cis"]
    if (cis_dt[,.N] > 1) {
      ld <- fread(sprintf("%s/%s.ld", cis_ld, aptamer))
      ld <- ld[, .(pos1=position1, pos2=position2, r2=correlation^2)]
      ld <- ld[pos1 %in% cis_dt[, pos]]
      ld <- ld[pos2 %in% cis_dt[, pos]]
      if (nrow(ld) > 1) { # 0 where LD between all variants < 0.001, so not stored in loaded ld table
				ld <- unique(rbind(ld, ld[,.(pos1=pos2, pos2=pos1, r2)]))
				ld[cis_dt, on = .(pos1=pos), ID1 := ALT_ID]
				ld[cis_dt, on = .(pos2=pos), ID2 := ALT_ID]
				ld <- as.matrix(dcast(ld, ID1 ~ ID2, value.var="r2", fill=0), rownames="ID1") 
				ld_mat[rownames(ld), colnames(ld)] <- ld
      }
    }
   
    # For trans-pQTLs we need to manually calculate LD
    trans_dt <- sub_dt[type == "trans"]
    if (trans_dt[, .N] > 1) {
      foreach(chrIdx = unique(trans_dt$chr), .errorhandling="remove") %do% {
        trans_chr_dt <- trans_dt[chr == chrIdx]
        if (trans_chr_dt[,.N] > 1) {
          # Use bgenix to extract each variant of interest
          foreach(rIdx = 1:nrow(trans_chr_dt)) %do% {
            cmd <- "bgenix"
            cmd <- paste(cmd, sprintf("-g data/INTERVAL/post_qc_data/imputed/uk10k_1000g_b37/imputed", chrIdx))
            cmd <- paste(cmd, sprintf("-incl-range %s:%s-%s", chrIdx, trans_chr_dt[rIdx, pos], trans_chr_dt[rIdx, pos]))
            cmd <- paste(cmd, sprintf("> %s/%s.bgen", ld_out, rIdx))
            system(cmd, wait=TRUE)
          }          

          # Combine bgen files 
          cmd <- sprintf("cat-bgen -g %s/*.bgen -og %s/combined.bgen", ld_out, ld_out)
          system(cmd, wait=TRUE)
      
					# Then use LDstore to calculate LD between variants
					cmd <- "ldstore"
					cmd <- paste(cmd, sprintf("--bcor %s/combined.bcor", ld_out))
					cmd <- paste(cmd, sprintf("--bgen %s/combined.bgen", ld_out))
					cmd <- paste(cmd, "--ld-thold 0 --n-threads 1")
					system(command=cmd, wait=TRUE)

					# Convert LDstore bcor binary into usable table
					cmd <- "ldstore"
					cmd <- paste(cmd, sprintf("--bcor %s/combined.bcor_1", ld_out))
					cmd <- paste(cmd, sprintf("--table %s/combined.ld", ld_out))
					system(command=cmd, wait=TRUE)
  
          # Read in table and process
          ld <- fread(sprintf("%s/combined.ld", ld_out))
					ld <- ld[, .(pos1=position1, pos2=position2, r2=correlation^2)]
					ld <- ld[pos1 %in% trans_chr_dt[, pos]]
					ld <- ld[pos2 %in% trans_chr_dt[, pos]]
					ld <- unique(rbind(ld, ld[,.(pos1=pos2, pos2=pos1, r2)]))
					ld[trans_chr_dt, on = .(pos1=pos), ID1 := ALT_ID]
					ld[trans_chr_dt, on = .(pos2=pos), ID2 := ALT_ID]
					ld <- as.matrix(dcast(ld, ID1 ~ ID2, value.var="r2"), rownames="ID1") 
					ld_mat[rownames(ld), colnames(ld)] <- ld
           
          # Clean up temporary files
          system(sprintf("rm %s/*.*", ld_out), wait=TRUE)
        }
      }
    }
    
    
    # Set diagonal to 1
    diag(ld_mat) <- 1

    # Write out
    write.table(ld_mat, sep=" ", quote=FALSE, file=sprintf("%s/instruments.ldmat", ld_out))
  }
}

