library(data.table)
library(foreach)
library(XML)

out_dir <- "analyses/pub/cardiometabolic_proteins/mouse_data"
dir.create(out_dir, showWarnings=FALSE)

# Load FDR significant associations
assocs <- fread("analyses/pub/cardiometabolic_proteins/all_assocs.tsv")
assocs <- assocs[Prot.FDR < 0.05]

# Split out to genes
genes <- assocs[, .(
  Gene=strsplit(Gene, "\\|")[[1]],
  EntrezID=as.integer(strsplit(Entrez.ID, "\\|")[[1]])),
  by=.(PRS, UniProt)]

# Load human to mouse gene mapping from MGI
MGI <- fread("data/HOM_MouseHumanSequence.rpt", check.names=TRUE)

# Get orthologs:
hg <- MGI[EntrezGene.ID %in% genes$EntrezID]
mm <- MGI[HomoloGene.ID %in% hg$HomoloGene.ID & Common.Organism.Name == "mouse, laboratory"]

hg <- hg[, .(HomoID = HomoloGene.ID, human.symbol = Symbol, EntrezID=EntrezGene.ID)]
mm <- mm[, .(HomoID = HomoloGene.ID, mouse.symbol = Symbol)]

orthologs <- merge(hg, mm, by="HomoID", all=TRUE)
genes <- merge(genes, orthologs, by = "EntrezID", all.x=TRUE)

# Validate genes with missing orthologs through the web interface of MGI
genes[Gene == "PDE4D", mouse.symbol := "Pde4d"] # http://www.informatics.jax.org/marker/MGI:99555
new_entry <- genes[Gene == "HBQ1"]
new_entry[, mouse.symbol := "Hbq1a"] # http://www.informatics.jax.org/marker/MGI:2685722
genes <- rbind(genes, new_entry)
genes <- unique(genes)

# Although genecards lists Pcdbh18 as the homolog of PCDBH10, MGI does not list PCDBH10 as 
# the human ortholog of Pcdbh18, and searching for PCDBH10 in MGI does not list any homologs
# (it returns Pcdbh10 which is a homolog for PCDBH5) and the NCBI gene page states PCDBH10 is
# conserved only in chimpanzees, rhesus monkey, and cow.

# The same applies for C5orf38/CEI and GGT2

# Manually validated cases where there is mouse symbol, but not transcript found in either HMDP or BxH.
# In these cases I used both MGI and NCBI gene to find and search for alternate symbols in the sysgen
# resource.
genes[, hmdp.symbol := mouse.symbol]
genes[, bxh.symbol := mouse.symbol]
genes[mouse.symbol == "Hbq1a", c("hmdp.symbol", "bxh.symbol") := "Hbq1"] # https://www.ncbi.nlm.nih.gov/gene/216635
genes[mouse.symbol == "C1s1",  c("hmdp.symbol", "bxh.symbol") := "C1s"] # https://www.ncbi.nlm.nih.gov/gene/50908
genes[mouse.symbol == "Ppp2r3a", bxh.symbol := "MMT00075066"] # sysgen resource lists transcripts, lookup by transcript id in gene info sheet
new_entry <- genes[mouse.symbol == "Ppp2r3a"]
new_entry[, bxh.symbol := "3222402P14Rik"] # alt symbol for Ppp2r3a with successful lookup in sysgen resource, matches genomic location
genes <- rbind(genes, new_entry)
genes[mouse.symbol == "Rida",  c("hmdp.symbol", "bxh.symbol") := "Hrsp12"] # alt symbol with successful lookup in sysgen resource
genes[mouse.symbol == "Try5", bxh.symbol := "1810049H19Rik"] # alt symbol 

# Obtain mapping of mouse genes to transcripts in each HMDP tissue:
hmdp_tissues <- c("HMDP_adipose", "HMDP_aorta", "HMDP_liver")
hmdp_transcripts <- foreach(tissue = hmdp_tissues, .combine=rbind) %do% {
  foreach(gene = unique(na.omit(genes$hmdp.symbol)), .combine=rbind) %do% {
    query <- paste0(
      "https://systems.genetics.ucla.edu/search",
      "?DATABASE_NAME=", tissue, 
      "&TABLE_NAME=GeneAnnotation",
      "&SEARCH_VARIABLE=gene_symbol",
      "&SEARCH_VALUE=", gene, 
      "&f=Submit")
    system(sprintf("wget -O %s/%s_%s_search.html '%s'", out_dir, tissue, gene, query), wait=TRUE)
    html <- htmlTreeParse(sprintf("%s/%s_%s_search.html", out_dir, tissue, gene))
    system(sprintf("rm %s/%s_%s_search.html", out_dir, tissue, gene), wait=TRUE)

    body <- xmlRoot(html)[["body"]]

    transcripts <- try({
      tbl <- body[["table"]]
      transcript_col <- NULL
      for (ii in seq_along(tbl[["thead"]][["tr"]])) {
        if (xmlValue(tbl[["thead"]][["tr"]][[ii]]) == "Transcript ID / Probeset ID") {
          transcript_col <- ii
        }
      }
      ts <- foreach(ii=2:length(tbl), .combine=c) %do% {
        xmlValue(tbl[[ii]][[transcript_col]])
      }
      ts
    }, silent=TRUE)

    if (class(transcripts) == "try-error") {
      transcripts <- NA_character_
    }
    data.table(gene, tissue, transcript=transcripts)
  }
}

# Iterate through transcripts tissue pairs and download available trait correlations
hmdp_cor <- foreach(rIdx = hmdp_transcripts[,.I], .combine=rbind) %do% {
  if (is.na(hmdp_transcripts[rIdx, transcript])) {
    return(NULL)
  } 
  tissue <- hmdp_transcripts[rIdx, tissue]
  gene <- hmdp_transcripts[rIdx, gene]
  transcript <- hmdp_transcripts[rIdx, transcript]
  query <- paste0(
    "https://systems.genetics.ucla.edu/correlation-ajax",
    "?dataset=", tissue,
    "&tid=", transcript,
    "&mode=trait", 
    "&cutoff=1.301", 
    "&f=Download")
  system(sprintf("wget -O %s/%s_%s_correlations.tsv '%s'", out_dir, tissue, transcript, query), wait=TRUE)
  trait_cor <- fread(sprintf("%s/%s_%s_correlations.tsv", out_dir, tissue, transcript))
  system(sprintf("rm %s/%s_%s_correlations.tsv", out_dir, tissue, transcript))
  trait_cor[, Tissue := tissue]
  trait_cor
}
setnames(hmdp_cor, c("N", "transcript", "symbol", "name", "phenotype", "corr", "pval", "tissue"))
hmdp_cor[hmdp_transcripts, on = .(tissue, transcript), gene := gene]
hmdp_cor <- hmdp_cor[, .(tissue, gene, transcript, phenotype, corr, pval)]

# Write out tables
fwrite(genes, sep="\t", quote=FALSE, file=sprintf("%s/mouse_orthologs.tsv", out_dir))
fwrite(hmdp_transcripts, sep="\t", quote=FALSE, file=sprintf("%s/hmdp_transcripts.tsv", out_dir))
fwrite(hmdp_cor, sep="\t", quote=FALSE, file=sprintf("%s/hmdp_trait_cor.tsv", out_dir))
