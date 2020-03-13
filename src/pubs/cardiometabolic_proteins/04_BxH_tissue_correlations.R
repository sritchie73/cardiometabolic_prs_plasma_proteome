library(data.table)
library(foreach)
library(WGCNA)

tissues <- c("adipose", "liver", "brain", "muscle")
genes <- fread("analyses/pub/cardiometabolic_proteins/mouse_data/mouse_orthologs.tsv", na.strings=c("", "NA"))

phen <- fread("data/BxH/cleaned/phenotype.txt")
phen <- melt(phen, id.vars=c("id", "sex"), variable.name="trait")
phen <- phen[!is.na(value)]

phen_info <- fread("data/BxH/cleaned/phenotype_info.txt")
gene_info <- fread("data/BxH/cleaned/expression_info.txt", colClasses=c("reporter_id"="character"))

assocs <- foreach(tn = tissues, .combine=rbind) %do% {
  expr <- fread(sprintf("data/BxH/cleaned/%s_expression.txt", tn), colClasses=c("array_id"="character"))
  expr <- melt(expr, id.vars="array_id", variable.name="sample_id")
  foreach(gn = na.omit(unique(genes$bxh.symbol)), .combine=rbind) %do% {
    ar_id <- gene_info[gene_sym == gn, reporter_id]
    foreach(an = ar_id, .combine=rbind) %do% {
      expr_gn <- expr[array_id == an]
      foreach(pn = phen_info$trait_name, .combine=rbind) %do% {
        phen_pn <- phen[trait == pn]
        dat <- merge(expr_gn, phen_pn, by.x="sample_id", by.y="id", suffixes=c(".expr", ".trait"))
        if (nrow(dat) < 2) return(NULL)
        dat[, .(tissue=tn, gene=gn, array_id=an, trait=pn,
								 Corr=bicorAndPvalue(value.expr, value.trait)[["bicor"]][1,1],
								 P=bicorAndPvalue(value.expr, value.trait)[["p"]][1,1],
								 N=bicorAndPvalue(value.expr, value.trait)[["nObs"]][1,1])]
      }
    }
  }
}
assocs[phen_info, on = .(trait=trait_name), trait := long_name]

fwrite(assocs, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/mouse_data/trait_associations.txt")
