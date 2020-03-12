library(data.table)
library(XML)
library(foreach)
library(doMC)
library(biomaRt)
library(openxlsx)

# Set up parallel environment
cpus <- as.numeric(system("echo $SLURM_CPUS_ON_NODE", intern=TRUE))
if (!is.na(cpus)) {
  registerDoMC(cpus)
} else {
  registerDoMC(1)
}

# make output directory
out_dir <- "analyses/processed_traits/somalogic_proteins/"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Load in post QC protein data
soma4k <- fread("data/INTERVAL/post_qc_data/phenotype/somalogic_proteomics/gwasqc/somalogic_qcgwas_4000.csv", colClasses=c("aliquot_id"="character"))

# Load in project 1074's omic identifier mapping sheet
id_map <- fread("data/INTERVAL/project_1074/omicsMap.csv", colClasses="character", na.strings=c("NA", ""))

# For this project, we'll use the genetic identifiers (IID) as the central sample identifiers.
soma4k[, join_key := aliquot_id]
soma4k <- id_map[,.(IID=Affymetrix_gwasQC_bl, soma4000_gwasQC_bl)][soma4k, on=.(soma4000_gwasQC_bl=join_key)]

# Split off the per-sample information
sample_info <- soma4k[, .(aliquot_id, soma4000_gwasQC_bl, IID, samplegroup, batch)]

# Keep only the IID column as sample identifiers
soma4k[, c("aliquot_id", "soma4000_gwasQC_bl", "samplegroup", "batch") := NULL]

# Drop samples that could not be mapped to IIDs
soma4k <- soma4k[!is.na(IID)]

# Transform to long format and write out
soma4k <- melt(soma4k, id.vars="IID")
fwrite(soma4k, file=sprintf("%s/traits.tsv", out_dir), sep="\t", quote=FALSE)

# Write out sample information and covariates file
fwrite(sample_info, file=sprintf("%s/sample_info.tsv", out_dir), sep="\t", quote=FALSE)
fwrite(sample_info[, .(IID, batch)], file=sprintf("%s/covariates.tsv", out_dir), sep="\t", quote=FALSE)

#-----------------------------------------------------------
# Load and update information sheet.
#-----------------------------------------------------------

# Load in information sheet
soma_info <- fread("data/INTERVAL/post_qc_data/phenotype/somalogic_proteomics/Documents/Somamer_info.tsv", na.strings=c("", "NA"))

# Load in updated information from SomaLogic. The new information is on their V4 platform,
# which (1) contains many new aptamers not included in the V3 release used in INTERVAL, 
# and importantly (2) only contains *1* representative aptamer per protein. In this release,
# they have de-duplicated aptamers, but will re-add these aptamers in the V5 release as they
# have apparently come full circle on the idea of duplicates (the thinking now is that these
# additional aptamers add information that may be information and that associations with 
# multiple aptamers for a given protein adds reliability to any association scan).
#
# In addition to the above, they have continued to test the specificity of each aptamer, and
# have found some of their original aptamers bind the fusion construct, not the protein itself.

# Read and excel sheet into a data.table
fread.xlsx <- function(...) {
  as.data.table(read.xlsx(...), check.names=TRUE)
}

v4_dir <- "data/INTERVAL/reference_files/SomaLogic_V4"

# Updated annotations in the V4 release. Contains many new aptamers not 
# measured in INTERVAL, and excludes a number of "duplicate aptamers" -
# those found by SomaLogic to have similar binding efficacy for the same 
# protein. Apparently these will be added back in in the next annotation
# release as they have re-evaluated duplicates (they add extra confirmatory
# information for association analyses).
v4 <- fread.xlsx(sprintf("%s/SOMAscan_Assay_v4_Annotations_version3.3.2.xlsx", v4_dir), sheet = 2)

# Reagents dropped from the V4 release as their aptamer was found to bind
# the "fusion construct" used to select the aptamer, not the human protein
# (double check terminology here) or aptamers "under investigation" - those
# still likely to bind a human protein.
V3notV4 <- fread.xlsx(sprintf("%s/V3 Delta V4 Reagents Update Aug 19 2019.xlsx", v4_dir))

# Filter V4 information to aptamers in the V3 sheet:
soma_update <- v4[SeqId %in% soma_info$Seq]

# Identify aptamers in the V3 sheet that are likely duplicate aptamers: those that
# are missing from the V4 sheet but likely measure the same target as one of the 
# aptamers in the V4 sheet. For matching by UniProt ID we need to make sure the format
# and sorting is consistent for aptamers with multiple UniProt IDs (protein complexes).
dups <- soma_info[!soma_update, on = .(Seq=SeqId)]
matched <- soma_info[Seq %in% soma_update$SeqId]
dups[, UniProtSorted := paste(sort(strsplit(UniProt, ",")[[1]]), collapse="|"), by=Seq]
soma_update[, UniProtSorted := paste(sort(strsplit(UniProt, "\\|")[[1]]), collapse="|"), by=SeqId]
v4[, UniProtSorted := paste(sort(strsplit(UniProt, "\\|")[[1]]), collapse="|"), by=SeqId]
matched[, UniProtSorted := paste(sort(strsplit(UniProt, ",")[[1]]), collapse="|"), by=Seq]
dups <- dups[UniProtSorted %in% v4$UniProtSorted | # Shares a UniProt ID (or IDs) with an apatamer in the V4 sheet
             UniProtSorted %in% matched$UniProtSorted | # Shares a UniProt ID (or IDs) with an aptamer in the V3 sheet 
             Target %in% v4$Target | # Shares a target name with an aptamer in the V4 sheet
             Target %in% matched$Target] # Shares a target name with an aptamer in the V3 sheet

# Annotate duplicates
dupby <- rbind(
  V4.UniProt = dups[UniProtSorted %in% v4$UniProtSorted, .(Seq, UniProtSorted, UniProt, Target)],
  V3.UniProt = dups[UniProtSorted %in% matched$UniProtSorted, .(Seq, UniProtSorted, UniProt, Target)],
  V4.Target = dups[Target %in% v4$Target, .(Seq, UniProtSorted, UniProt, Target)],
  V3.Target = dups[Target %in% matched$Target, .(Seq, UniProtSorted, UniProt, Target)],
  idcol="dupby")

# The following table is empty: this means there are no aptamers for associated 
# with > 1 Target if we were to take the information from matches by UniProt or
# Target name in both V3 and V4 annotations.
unique(dupby, by=c("Seq", "UniProt", "Target"))[,.N,by="Seq"][N > 1]

# Organise dupby annotations to take 1 unique entry per aptamer. When merging
# duplicates into the updated table we may need to handle those matched by V3 
# and those matched by V4 information differently
dupby[, dupby := factor(dupby, levels=c("V4.Target", "V4.UniProt", "V3.Target", "V3.UniProt"))]
dupby <- dupby[order(dupby)]
dupby <- dupby[, .SD[1],  by=Seq] # i.e. best match per aptamer

# Each dupby case needs to be handled differently when adding to the updated annotation sheet,
# but in each case we must make sure that aptamer specific information is not copied across.
tomerge <- dupby[dupby == "V4.Target", .(Seq, Target)]
tomerge <- tomerge[v4, on = .(Target), nomatch=0]
tomerge[, c("Seq", "SeqId") := .(NULL, Seq)] # SeqId column should have Seq from dupby.
tomerge[, c("Type", "Characterization.Info", "Mass.Spec.Confirmation.in.Matrix",
            "Total.CV.Plasma", "Total.CV.Serum") := NULL] # Aptamer specific information
soma_update <- rbind(soma_update, tomerge, fill=TRUE)

tomerge <- dupby[dupby == "V4.UniProt", .(Seq, UniProtSorted)] 
tomerge <- tomerge[v4, on = .(UniProtSorted), nomatch=0]
tomerge[, c("Seq", "SeqId") := .(NULL, Seq)] # SeqId column should have Seq from dupby.
tomerge[, c("Type", "Characterization.Info", "Mass.Spec.Confirmation.in.Matrix",
            "Total.CV.Plasma", "Total.CV.Serum") := NULL] # Aptamer specific information
soma_update <- rbind(soma_update, tomerge, fill=TRUE)

tomerge <- dupby[dupby == "V3.Target", .(Seq, Target)]
tomerge <- tomerge[matched, on = .(Target), .(Seq, V4.SeqMatch = i.Seq), nomatch=0]
tomerge <- tomerge[soma_update, on = .(V4.SeqMatch=SeqId), nomatch=0]
tomerge[, c("SeqId", "Seq", "V4.SeqMatch") := .(Seq, NULL, NULL)] # SeqId column should have Seq from dupby.
tomerge[, c("Type", "Characterization.Info", "Mass.Spec.Confirmation.in.Matrix",
                        "Total.CV.Plasma", "Total.CV.Serum") := NULL] # Aptamer specific information
soma_update <- rbind(soma_update, tomerge, fill=TRUE)

# Add annotation to each aptamer indicating match type
soma_update[, update_anno := "V3_in_V4"]
soma_update[is.na(Type), update_anno := "V3_duplicate_in_V4"]

# There are some duplicate targets in the V4 sheet for these duplicate aptamers. In all
# but one case these are because there are multiple aptamers in the V4 sheet with the 
# same target but slightly different entries in the Apparent.Kd.(M) column (Estimated 
# apparent Kd (molar units) based on 4pl fit to buffer standard curve with the target 
# protein(s) used in the selection of the SOMAmer reagent). In these cases we just set
# the value of this column to NA.
multi <- soma_update[update_anno == "V3_duplicate_in_V4",.N,by=SeqId][N>1]
soma_update[SeqId %in% multi$SeqId, `Apparent.Kd.(M)` := NA]
soma_update <- unique(soma_update)

# Two entry remains - in each case it is because there are two aptamers with the same
# Target, but one is missing information from UniProt. Therefore, we just take the entry
# with the UniProt information.
multi <- soma_update[update_anno == "V3_duplicate_in_V4",.N,by=SeqId][N>1]
soma_update <- soma_update[!(SeqId %in% multi$SeqId & UniProt == "")]

# No longer need UniProtSorted column
soma_update[, UniProtSorted := NULL]

# Identify aptamers in the V3 sheet that are missing from the V4 sheet, and could not
# be matched by UniProt or Target to another aptamer in included in the V3 or V4 sheet.
v3_missing <- soma_info[!soma_update, on = .(Seq=SeqId)]

# Add these, using the information present in the V3 Sheet
tomerge <- v3_missing[, .(SeqId=Seq, Target, TargetFullName, UniProt, Organism, 
                          Gene.Name.Name=GENE, update_anno = "V3_missing_in_V4")]
soma_update <- rbind(soma_update, tomerge, fill=TRUE)

# Add annotation information sent by somalogic about V3 aptamers not included in the V4 sheet:
soma_update <- merge(soma_update, V3notV4, by.x = "SeqId", by.y="SeqID.in.v3.i5K", all.x=TRUE)
soma_update[!is.na(Analyte) & Annotation == "MOUSE_Fc", Target := "Fc_MOUSE"]
soma_update[!is.na(Analyte) & Annotation == "MOUSE_Fc",
            c("TargetFullName", "UniProt", "Organism", "Type", "Characterization.Info",
              "Mass.Spec.Confirmation.in.Matrix", "UniProt.Id.Current.at.Uniprot",
              "UniProt.Full.Name", "Uniprot.Short.Names", "Uniprot.Alternate.Names",
              "Reactome.Pathways", "Subcellular.Locations", "Gene.Name.Name", 
              "Ensembl.Gene.ID", "Entrez.Gene.ID", "HGNC.ID", "GO.Function.Annotations", 
              "GO.Process.Annotations", "GO.Component.Annotations") := NA]
soma_update[!is.na(Analyte) & Target != "Fc_MOUSE" & Type != "Deprecated", Type := Annotation]
soma_update[, c("Analyte", "Annotation") := NULL]

# We were also sent a sheet containing information about pQTLs changed, but it contains
# no additional information useful here.

# INTERVAL specific information to add back in
soma_info <- soma_info[, .(Seq, SOMAMER_ID, Plex_1129, 
                           INTERVAL_Batch1_SOMALOGIC_QC=Batch1_SOMALOGIC_QC,
                           INTERVAL_Batch2_SOMALOGIC_QC=Batch2_SOMALOGIC_QC,
                           INTERVAL_Batch1_CV=Batch1_CV,
                           INTERVAL_Batch2_CV=Batch2_CV)]

# Add identifiers that match the column names in the 4K sheet. In most cases
# these are constructed via the gsub statement below, but a handful need 
# correcting. Note that column order in:
# data/INTERVAL/post_qc_data/phenotype/somalogic_proteomics/gwasqc/somalogic_qcgwas_4000.csv
# is the same as the row order in:
# data/INTERVAL/post_qc_data/phenotype/somalogic_proteomics/Documents/Somamer_info.tsv
# which allows us to match the four identifiers needing correcting.
soma_info[,variable := gsub("\\.", "", tolower(SOMAMER_ID))] 
soma_info[variable == "ighg1ighg2ighg3ighg4igkigl2744572", variable := "ighg1ighg2ighg3ighg4igkigl274457"]
soma_info[variable == "ighg1ighg2ighg3ighg4igkigl3700154", variable := "ighg1ighg2ighg3ighg4igkigl370015"]
soma_info[variable == "ywhabywhaeywhagywhahywhaqywhazsfn4179573", variable := "ywhabywhaeywhagywhahywhaqywhazsf"]
soma_info[variable == "ywhabywhaeywhagywhahywhaqywhazsfn4707502", variable := "v2217"]

# Add into the updated somamer information sheet
soma_update <- merge(soma_info, soma_update, by.x="Seq", by.y="SeqId")
setnames(soma_update, "Seq", "SeqId")

# The V4 sheet contains the updated full, short, and alternate names listed on UniProt. For
# aptamers listed only on the V3 sheet, we need to pull down this information.

# Pull down from UniProt the recommended names and primary gene symbol associated
# with each UniProtKB identifer for aptamers in the V3 sheet but not the V4 sheet,
# which have up to date information.
uids <- soma_update[update_anno == "V3_missing_in_V4" & Organism == "Human", UniProt]
uids <- unlist(strsplit(uids, ","))
uids <- unique(uids)
uids <- na.omit(uids)
uniprot_info <- foreach(uid = uids, .combine=rbind, .errorhandling="remove") %dopar% {
  system(sprintf("wget -P %s https://www.uniprot.org/uniprot/%s.xml", out_dir, uid), wait=TRUE)
  xml <- xmlTreeParse(sprintf("%s/%s.xml", out_dir, uid))
  system(sprintf("rm %s/%s.xml", out_dir, uid), wait=TRUE)
  
  root <- xmlRoot(xml)

  # Get the entry name and strip out the "_HUMAN"
  entry_name <- xmlValue(root[["entry"]][["name"]])
  entry_name <- gsub("_HUMAN", "", entry_name)
  
  # Get the recommended full name of the protein, taking the first entry if multiple 
  # (unclear if this is ever the case)
  rn_node <- root[["entry"]][["protein"]][["recommendedName"]]
  full_name <- xmlValue(rn_node[["fullName"]])

  # Get the recommended short name. If there are multiple, choose the shortest one.
  sn_nodes <- xmlElementsByTagName(rn_node, "shortName")
  short_names <- foreach(node = sn_nodes, .combine=c) %do% {
    xmlValue(node)
  }
  if (is.null(short_names)) {
    short_name <- NA_character_
  } else {
    short_name <- short_names[which.min(nchar(short_names))]
  }

  # Get the primary gene symbol
  gn_nodes <- root[["entry"]][["gene"]]
  gene <- foreach(node = gn_nodes, .combine=c) %do% {
    if (xmlAttrs(node)["type"] == "primary") {
      return(xmlValue(node))
    }
  }

  data.table(UniProt=uid, EntryName=entry_name, RecommendedFull=full_name, 
             RecommendedShort=short_name, SYMBOL=gene)
}

soma_update[uniprot_info, on = .(UniProt), 
            c("UniProt.Id.Current.at.Uniprot", "UniProt.Full.Name", "Uniprot.Short.Names", "Gene.Name.Name") := 
            .(UniProt, RecommendedFull, RecommendedShort, SYMBOL)]

# To obtain the genomic location we will look up each gene by Entrez ID via
# the NCBI gene website. Some entries are missing Entrez IDs, primarily 
# aptamers present in the V3 sheet but missing from the V4 sheet.

# Lookup Entrez IDs for each UniProt Identifier using biomart
missing_eid <- soma_update[(is.na(Entrez.Gene.ID) | Entrez.Gene.ID == "")]
missing_eid <- missing_eid[Organism == "human" & Target != "Fc_MOUSE" & UniProt != "" & !is.na(UniProt)]
missing_eid <- missing_eid[,.(SeqId, UniProt)]
missing_eid <- missing_eid[,.(UniProt=strsplit(UniProt, "\\|")[[1]]), by=SeqId]

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
bm <- as.data.table(getBM(mart = ensembl, values = missing_eid$UniProt, 
            filters="uniprot_gn_id", attributes=c("uniprot_gn_id", "entrezgene_id")))
setnames(bm, c("UniProt", "ENTREZ"))
bm[, ENTREZ := as.character(ENTREZ)]

missing_eid[bm, on = .(UniProt), Entrez.Id := ENTREZ]

# 14 UniProt IDs did not have an associated Entrez ID in biomart. We enter these manually by
# looking up their gene symbol on the UniProt page, then searching NCBI gene for the gene:
missing_eid[is.na(Entrez.Id), Entrez.Id := c("6927", "23269", "5078", "29984",
                                             "7004", "6000", "55786", "7005",
                                             "3501", "3503", "51375", "9674")] 
missing_eid <- missing_eid[order(UniProt)][order(SeqId)]
missing_eid <- missing_eid[, .(UniProt=paste(UniProt, collapse="|"), 
                               Entrez.Id=paste(Entrez.Id, collapse="|")), 
                           by=SeqId]

soma_update[missing_eid, on = .(SeqId), Entrez.Gene.ID := Entrez.Id]

# I extract the genomic location for each Entrez ID directly from the webpage on NCBI
# gene because unlike other databases, it provides me exactly one genomic location.
eids <- soma_update[!is.na(Entrez.Gene.ID) & Entrez.Gene.ID != "", Entrez.Gene.ID]
eids <- unique(unlist(strsplit(eids, "\\|")))
ncbi_info <- foreach(eid = eids, .combine=rbind, .errorhandling="remove") %dopar% {
  system(sprintf("wget -P %s https://www.ncbi.nlm.nih.gov/gene/%s", out_dir, eid), wait=TRUE)
  xml <- xmlTreeParse(sprintf("%s/%s", out_dir, eid))
  system(sprintf("rm %s/%s", out_dir, eid), wait=TRUE)

  root <- xmlRoot(xml)

  # Traverse tree to get to node corresponding to content of the web page. Inspect element in
  # chrome was used to guide node number selection at each step.
  node <- root[[2]][[1]][[1]][[1]][[1]]  # <body ...> <div ...> <div ...> <form ...> <div ...>
  node <- node[[7]] # <div id="maincontent">
  node <- node[[1]][[8]] # <div ...> <div id="rprt full-rprt">
  content <- node[[2]] # <div class="rprt-body">

  # Extract the recommended gene symbol provided by HGNC
  node <- content[[1]] # <div class="rprt-section gene-summary">
  node <- node[[2]][[1]][[1]] # <div ...> <div ...> <div ...>
  symbol <- xmlValue(node[[2]][[1]])

  # Extract HTML node corresponding to the genomic location table
  node <- content[[2]] # <div class="rprt-section gene-genomic-context">
  node <- node[[2]] # <div class="rprt-section-body...>
  tbl <- node[[1]][[3]][[1]][[1]] # <div ...> <div class="uincbigrid-outer-div"> <div class="uincbigrid-inner-div"> <table ...>

  # Parse HTML table to data.table
  header <- sapply(tbl[[1]][[1]], xmlValue)
  rows <- foreach(rIdx = seq_along(tbl[[2]]), .combine=rbind) %do% {
    foreach(cIdx = seq_along(tbl[[2]][[rIdx]]), .combine=cbind) %do% {
      xmlValue(tbl[[2]][[rIdx]][[cIdx]])
    }
  }
  tbl <- as.data.table(rows)
  setnames(tbl, header)

  # Extract genomic location:

  # In some instances more rows than expected are found. These occur when
  # (1) a gene has a position on both the X and Y chromosomes
  # (2) the gene only occurs on alternative haplotypes to the primary
  #     reference genome (e.g. for some HLA genes).
  # In the first case, we just take the position on the X chromosome,
  # in the second we take the smallest alternative haplotype (by number).
  if (length(unique(tbl$Chr)) > 1) {
    if (all(sort(unique(tbl$Chr)) == c("X", "Y"))) {
      tbl <- tbl[Chr == "X"]
    } else {
      altnum <- tbl[, Chr]
      chr <- gsub(" .*", "", altnum)
      altnum <- gsub(".*ALT_REF_LOCI_", "", altnum)
      altnum <- gsub(")$", "", altnum)
      minalt <- min(altnum)
      tbl <- tbl[Chr == sprintf("%s (ALT_REF_LOCI_%s)", chr, minalt)]
      tbl[, Chr := chr]
    }
  }

  if (any(grepl("GRCh37", tbl$Assembly))) {
    row <- tbl[Assembly %like% "GRCh37"]
    build <- "GRCh37"
  } else {
    # Sometimes a position is not available on GRCh37, i.e.
    # where the gene's position was unknown at the time of that build.
    row <- tbl[Assembly %like% "GRCh38"]
    build <- "GRCh38"
  }

  chr <- row[, Chr]
  locstring <- row[, Location]

  start <- gsub(".*\\(", "", locstring)
  start <- gsub("\\.\\..*", "", start)
  end <- gsub(".*\\.\\.", "", locstring)
  end <- gsub("[^0-9]*$", "", end)

  strand <- sapply(locstring, function(x) {
    if (grepl("[0-9])$", x)) {
      "+"
    } else if (grepl("complement", x)) {
      "-"
    } else {
      "*"
    }
  })

  data.table(ENTREZ = eid, SYMBOL=symbol, chr, start, end, build, strand)
}

# Link information back to main table
eidmap <- soma_update[!is.na(Entrez.Gene.ID) & Entrez.Gene.ID != "", 
                      .(EID = strsplit(Entrez.Gene.ID, "\\|")[[1]]), 
                      by=.(SeqId, Entrez.Gene.ID)]
ncbi_info <- ncbi_info[eidmap, on = .(ENTREZ = EID)]

# There are 23 EIDs where information is missing or the lookup failed
# in some way, we will fill these out manually:
ncbi_info[ENTREZ == "3493", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHA1", "14", "106173505", "106175001", "GRCh37", "-")]
ncbi_info[ENTREZ == "3494", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHA2", "14", "106053274", "106054731", "GRCh37", "-")]
ncbi_info[ENTREZ == "3495", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHD", "14", "106304737", "106312010", "GRCh37", "-")]
ncbi_info[ENTREZ == "3497", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHE", "14", "106066403", "106068064", "GRCh37", "-")]
ncbi_info[ENTREZ == "3500", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHG1", "14", "106207810", "106209407", "GRCh37", "-")]
ncbi_info[ENTREZ == "3501", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHG2", "14", "106109540", "106111126", "GRCh37", "-")]
ncbi_info[ENTREZ == "3502", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHG3", "14", "106232251", "106237742", "GRCh37", "-")]
ncbi_info[ENTREZ == "3503", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHG4", "14", "106090813", "106092402", "GRCh37", "-")]
ncbi_info[ENTREZ == "3507", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGHM", "14", "106318298", "106322322", "GRCh37", "-")]
ncbi_info[ENTREZ == "3535", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGL", "22", "22380474", "23265085", "GRCh37", "+")]
ncbi_info[ENTREZ == "3803", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("KIR2DL2", "19", "21890", "36439", "GRCh37", "+")]
ncbi_info[ENTREZ == "3813", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("KIR2DL2", "19", "21890", "36439", "GRCh37", "+")]
ncbi_info[ENTREZ == "4276", # Replaced by entry 100507436
          c("ENTREZ", "SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("100507436", "MICA", "6", "31367561", "31383090", "GRCh37", "+")]
ncbi_info[ENTREZ == "5645", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("PRSS2", "7 (PATCHES)", "974016", "977604", "GRCh37", "+")]
ncbi_info[ENTREZ == "11026", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("LILRA3", "19", "54799854", "54804625", "GRCh37", "-")]
# This one links to a mouse gene, Casq1. Looking up soma_update for entries with
# this Entrez.Id, we see two proteins and genes. One for CGA, which has the correct
# Entrez.Id, and one for TSHB, which does not have an entrez Id in this field. We
# therefore overwrite this mouse gene entry with the NCBI gene entry for TSHB
ncbi_info[ENTREZ == "12372",
          c("ENTREZ", "SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("7252", "TSHB", "1", "115572445", "115576930", "GRCh37", "+")]
# Multiple positions listed, use ; as a field separator to distinguish from | when
# there are multiple genes for an apatamer
ncbi_info[ENTREZ == "50802", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("IGK", "2", "89890568;89156874", "90274235;89630436", "GRCh37", "+;-")]
ncbi_info[ENTREZ == "57292", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("KIR2DL5A", "19", "86691", "96155", "GRCh37", "+")]
# NCBI search linked to right record, but page didnt parse. Listed as not in annotation
# release, with no genomic coordinates beyond 21p11.2.
# UCSC genome browser and RefSeq beg to differ:
# BAGE3 at chr21:11020842-11098925 - (NM_182481) B melanoma antigen 3 precursor
ncbi_info[ENTREZ == "85318", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("BAGE3", "21", "11020842", "11098925", "GRCh37", "-")]
ncbi_info[ENTREZ == "388372", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("CCL4L1", "17", "34640034", "34641841", "GRCh37", "+")]
ncbi_info[ENTREZ == "652493", SYMBOL := "LOC652493"] # Status: WITHDRAWN 
ncbi_info[ENTREZ == "1724716", SYMBOL := "HIV2gp6"] # Non-human target
ncbi_info[ENTREZ == "100132285", 
          c("SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("KIR2DS2", "19", "131449", "145745", "GRCh37", "+")]

# Validate entries with non-standard chromosomes using UCSC genome browser/refseq:
ncbi_info[SYMBOL == "PRSS2", # UCSC genome browser under alternate symbol TRYP2
          c("chr", "start", "end", "build", "strand") := 
          .("7", "141968101", "141972068", "GRCh37", "-")]
ncbi_info[ENTREZ == "107987462", # soma_update lists LILRB3 as gene, fix
          c("ENTREZ", "SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("11025", "LILRB3", "19", "54720147", "54726997", "GRCh37", "-")]
ncbi_info[SYMBOL == "CRYAA2", # UCSC genome browser under alternate symbol CRYAA
          c("chr", "start", "end", "build", "strand") := 
          .("21", "44589141", "44592913", "GRCh37", "+")]
# RAB7B not aligned to GRCh37
ncbi_info[ENTREZ == "107080644", # soma_update lists GNMT as gene, fix
          c("ENTREZ", "SYMBOL", "chr", "start", "end", "build", "strand") :=
          .("27232", "GNMT", "6", "42928500", "42931618", "GRCh37", "+")]

# Note the right join above ensures the row order of EIDs is the same
# as the order in Entrez.Gene.ID, so when we collapse the information
# the order of genes is preserved across columns (when there are 
# multiple genes separated by a | ).
ncbi_info <- ncbi_info[,.(
  Entrez.Gene.ID = paste(ENTREZ, collapse="|"),
  Gene.Name = paste(SYMBOL, collapse="|"),
  chr = paste(chr, collapse="|"),
  start = paste(start, collapse="|"),
  end = paste(end, collapse="|"),
  build = ifelse(length(unique(build)) == 1, unique(build), paste(build, collapse="|")),
  strand = paste(strand, collapse="|")),
  by = SeqId]


# Merge info:
soma_update[, Entrez.Gene.ID := NULL] # some have been updated manually due to verification checks
soma_update <- merge(soma_update, ncbi_info, by="SeqId", all.x=TRUE)

# There are 149 entries where the new gene symbol does not match the one provided by somalogic.
# Manually verify and fix these where necessary, using the UniProt ID as ground truth for target.
# In the vast majority of cases, these are just out of date / alternate gene symbols provided 
# by SomaLogic or UniProt.

# Entrez ID given by somalogic does not match gene listed on UniProt page.
soma_update[SeqId == "10600-24", 
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("1047", "CLGN", "4", "141309607", "141348815", "GRCh37", "-")]
soma_update[SeqId == "11273-176",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("653689", "GSTT2B", "22", "24299601", "24303416", "GRCh37", "-")]
soma_update[SeqId == "11626-7",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("7335", "UBE2V1", "20", "48697661", "48732496", "GRCh37", "-")]
soma_update[SeqId == "13944-3",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("6818", "SULT1A3", "16", "30210549", "30215650", "GRCh37", "+")]
soma_update[SeqId %in% c("14623-26", "9169-14"),
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("6612", "SUMO3", "21", "46225532", "46238044", "GRCh37", "-")]
soma_update[SeqId == "5097-14",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("3813", "KIR3DS1", "19", "70097", "84658", "GRCh37", "+")]
soma_update[SeqId == "5726-49",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("653423", "SPAG11A", "8", "7705209", "7721319", "GRCh37", "+")]
soma_update[SeqId == "7059-14",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("79168", "LILRA6", "19", "54740469", "54746724", "GRCh37", "-")]
soma_update[SeqId == "7082-2",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("192134", "B3GNT6", "11", "76745385", "76753020", "GRCh37", "+")]
soma_update[SeqId == "7895-108",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("56246", "MRAP", "21", "33664124", "76753020", "GRCh37", "+")]
soma_update[SeqId == "7982-10",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("64983", "MRPL32", "7", "42971804", "42977456", "GRCh37", "+")]
soma_update[SeqId == "9264-11",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("1519", "CTSO", "4", "156845270", "156875048", "GRCh37", "-")]
soma_update[SeqId == "9867-23",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("8789", "FBP2", "9", "97321002", "97356075", "GRCh37", "-")]

# Multiple genes corresponding to single uniprot
soma_update[SeqId == "2783-18",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("6349|414062", "CCL3L1|CCL3L3", "17|17", "34623842|34522268", 
              "34625730|34524147", "GRCh37", "-|-")]
soma_update[SeqId == "7842-52",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("8293|728492", "SERF1A|SERF1B", "5|5", "70196490|69321078", 
              "70214357|69338935", "GRCh37", "+|+")]

# Update genes based on cross-reference checking to the extensive manual checks I performed
# when updating the V3 platform:

# P04745 is encoded by AMY1A, AMY1B, and AMY1C not just AMY1A:
soma_update[UniProt == "P04745",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("276|277|278", "AMY1A|AMY1B|AMY1C", "1|1|1", "104193695|104230039|104287833",
              "104207174|104243519|104301312", "GRCh37", "+|-|+")]
# Q9Y6F8 is encoded by CDY1 and CDY1B not just CDY1:
soma_update[UniProt == "Q9Y6F8",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("9085|253175", "CDY1|CDY1B", "Y|Y", "27768264|26191376", "27771049|26194161",
              "GRCh37", "+|-")]
# SomaLogic lists P01233 as a UniProt identifier for hCG, but this is obselete, replaced by
# both P0DN86 and P0DN87. Only P0DN86 is listed as a current uniprot identifier for this
# aptamer, but our V3 sheet lists P0DN87 as well, so we will add this and all gene information.
soma_update[UniProt == "P01215|P01233",
            c("UniProt.Id.Current.at.Uniprot", "Entrez.Gene.ID", "Gene.Name", 
              "chr", "start", "end", "build", "strand") :=
            .("P01215|P0DN86|P0DN87", "1081|1082|93659|94027|94115", "CGA|CGB3|CGB5|CGB7|CGB8",
              "6|19|19|19|19", "87795216|49526126|49547102|49557531|49550895",
              "87804865|49527632|49548568|49561603|49552368", "GRCh37", "-|-|+|-|-")]
# CRYAA is the correct entry for P02489, not CRYAA2
soma_update[UniProt == "P02489",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("1409", "CRYAA", "21", "44589141", "44592915", "GRCh37", "+")]
# CPN2 is the correct entry for P22792, not CYP11B2
soma_update[UniProt == "P22792",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("1370", "CPN2", "3", "194060494", "194072063", "GRCh37", "-")]
# P59665 is coded by DEFA1 and DEFA1B
soma_update[UniProt == "P59665",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("1667|728358", "DEFA1|DEFA1B", "8|8", "6835171|6854288", "6837614|6856724",
              "GRCh37", "-|-")]
# Q8WTQ1 is coded by DEFB104A and DEFB104B
soma_update[UniProt == "Q8WTQ1",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("140596|503618", "DEFB104A|DEFB104B", "8|8", "7693993|7327830", "7698764|7332604",
              "GRCh37", "+|-")]
# Q8IZN7 is coded by DEFB107A and DEFB107B
soma_update[UniProt == "Q8IZN7",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("245910|503614", "DEFB107A|DEFB107B", "8|8", "7669242|7353368", "7673238|7366833",
              "GRCh37", "-|+")]
# P30042 listed by SomaLogic is obselete, replaced by A0A0B4J2D5 and P0DPI2. Only the former
# is listed as the current UniProt Identifier, but the Entrez ID provided by UniProt links to
# GATD3A, which corresponds to P0DPI2 rather than GATD3B which corresponds to A0A0B4J2D5. The
# entry for GATD3A also corresponds to the Target Full Name provided, so we modify the UniProt
# information to match this target entry
soma_update[UniProt == "P30042", 
            c("UniProt.Id.Current.at.Uniprot", "UniProt.Full.Name", "Uniprot.Short.Names", "Uniprot.Alternate.Names") :=
            .("P0DPI2", "Glutamine amidotransferase-like class 1 domain-containing protein 3A, mitochondrial", NA, NA)]
# P68431 is encoded by ten genes:
soma_update[UniProt == "P68431",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("8350|8358|8352|8351|8353|8968|8355|8357|8354|8356",
              "H3C1|H3C2|H3C3|H3C4|H3C6|H3C7|H3C8|H3C10|H3C11|H3C12",
              "6|6|6|6|6|6|6|6|6|6",
              "26020718|26031817|26045639|26196982|26224427|26250370|26271146|27777842|27839623|27858093",
              "26021186|26032288|26046097|26199521|26227706|26250835|26271612|27778314|27840099|27858570",
              "GRCh37", "+|-|+|-|+|-|-|+|-|-")]
# P47929 is coded by LGALS7 and LGALS7B
soma_update[UniProt == "P47929",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("3963|653499", "LGALS7|LGALS7B", "19|19", "39261608|39279850", "39264157|39282394",
              "GRCh37", "-|+")]
# Q9NQU5 is coded by PAK6 and BUB1B-PAK6 readthrough
soma_update[UniProt == "Q9NQU5",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("56924|106821730", "PAK6|BUB1B-PAK6", "15|15", "40509629|40217428", "40569688|40277487",
              "GRCh37", "+|+")]
# Q96QH8 is coded by SPACA5 and SPACA5B
soma_update[UniProt == "Q96QH8",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("389852729201", "SPACA5|SPACA5B", "X|X", "47863734|47990039", "47869130|47991995",
              "GRCh37", "+|+")]

# Drop old gene.name column
soma_update[, Gene.Name.Name := NULL]

# Add in cross-reactivity testing from the pQTL paper for aptamers not in the V4 release.
v3_cr <- fread.xlsx("data/INTERVAL/reference_files/SunB_2018_pQTL_study_supp_tables.xlsx", sheet=3, startRow=3)
v3_cr[, SeqId := gsub("_.*", "", SOMAmer.ID)]
v3_cr <- v3_cr[SeqId %in% soma_update[update_anno != "V3_in_V4", SeqId]]

v3_cr <- v3_cr[, .(Tested=paste0(paste(Related.Protein, collapse=", "), ".")), by=.(SeqId, Binding.result)]
v3_cr[, Tested := paste(Binding.result, "for", Tested)]
v3_cr <- v3_cr[, .(CrossReact = paste(Tested, collapse=" ")), by=SeqId]

soma_update[v3_cr, on = .(SeqId), Characterization.Info := CrossReact] 

# Add or update annotations in Type column for filtering:
soma_update[Target %in% c("Tenascin", "HSP 70"), Type := "Contaminant"]
soma_update[Target == "Fc_MOUSE", Type := "Fusion Construct"]
soma_update[Organism == "Human", Organism := "human"]
soma_update[(is.na(Type) | Type == "Protein") & Organism != "human", Type := "Non-Human"]
soma_update[is.na(Type), Type := "Protein"]

# Add gene information for remaining UniProt identifiers with blank Entrez.Gene.ID:
soma_update[UniProt == "P06753",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("7170", "TPM3", "1", "154127780", "154164611", "GRCh37", "-")]  
soma_update[UniProt == "P06732",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("1158", "CKM", "19", "45809671", "45826302", "GRCh37", "-")]
soma_update[UniProt == "P04150",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("2908", "NR3C1", "5", "142657496", "143113322", "GRCh37", "-")]
soma_update[UniProt == "Q92794",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("7994", "KAT6A", "8", "41786997", "41909506", "GRCh37", "-")]
soma_update[UniProt == "Q05655",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("5580", "PRKCD", "3", "53195223", "53226733", "GRCh37", "+")]
soma_update[UniProt == "P03372",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("2099", "ESR1", "6", "152011631", "152424409", "GRCh37", "+")]
soma_update[UniProt == "O75582",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("9252", "RPS6KA5", "14", "91335086", "91526993", "GRCh37", "-")]
soma_update[UniProt == "P47710",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("1446", "CSN1S1", "4", "70796799", "70812288", "GRCh37", "-")]
soma_update[UniProt == "P36894",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("657", "BMPR1A", "10", "88516396", "88684945", "GRCh37", "+")]
soma_update[UniProt == "Q93038",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("8718", "TNFRSF25", "1", "6520787", "6527432", "GRCh37", "-")]
soma_update[UniProt == "Q9NQC3",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("57142", "RTN4", "2", "55199325", "55307770", "GRCh37", "-")]
soma_update[UniProt == "P01116",
            c("Entrez.Gene.ID", "Gene.Name", "chr", "start", "end", "build", "strand") :=
            .("3845", "KRAS", "12", "25358180", "25403870", "GRCh37", "-")]

# Write out protein information
fwrite(soma_update, file=sprintf("%s/trait_info.tsv", out_dir), sep="\t", quote=FALSE)

# Note the user should prefer the SomaLogic name because there are many
# instances where the its targetting an isoform of the protein - but 
# each protein/aptamer of interest should be validated to make sure the
# name is not just an out of date UniProt name.
