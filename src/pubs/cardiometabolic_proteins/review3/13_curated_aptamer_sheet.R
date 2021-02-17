library(data.table)
library(openxlsx)

# Load Aptamer information
sinfo <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")
sinfo <- sinfo[, .(variable, SOMAMER_ID, Aptamer=SeqId, Type, IntendedTarget=TargetFullName, 
                   UniProt=UniProt.Id.Current.at.Uniprot, Gene=Gene.Name, 
                   Entrez=Entrez.Gene.ID, Chr=chr, Start=start, End=end, Cross_Reactivity=Characterization.Info, 
                   Mass_Spec_Confirmation=Mass.Spec.Confirmation.in.Matrix)]

# Fix bad entries (Aptamers for the same target with different/missing gene/uniprot information)
sinfo[IntendedTarget == "14-3-3 protein family", UniProt := "P61981|Q04917"]
sinfo[IntendedTarget == "Induced myeloid leukemia cell differentiation protein Mcl-1",
  c("UniProt", "Gene", "Entrez", "Chr", "Start", "End") :=
  .("Q07820", "MCL1", "4170", "1", "150547027", "150552214")]
sinfo[IntendedTarget == "Protein delta homolog 1",
  c("UniProt", "Gene", "Entrez", "Chr", "Start", "End") :=
  .("P80370", "DLK1", "8788", "14", "101193202", "101201467")]
sinfo[IntendedTarget == "Stromal cell-derived factor 1",
  c("UniProt", "Gene", "Entrez", "Chr", "Start", "End") :=
  .("P48061", "CXCL12", "6387", "10", "44865601", "44880545")]

# Curate aptamer sensitivity and specificity information
emilsson <- rbind(fill=TRUE,
  as.data.table(read.xlsx("data/Emilsson_etal_2018/supp_tables.xlsx", sheet="Table S3", startRow=2)),
  as.data.table(read.xlsx("data/Emilsson_etal_2018/supp_tables.xlsx", sheet="Table S4", startRow=2)))
soma_by_gene <- sinfo[,.(Gene=strsplit(Gene, "\\|")[[1]]),by=Aptamer]
soma_by_gene[emilsson, on = .(Gene=Gene.Symbol), Mass_Spec_Confirmation := TRUE]
soma_by_gene <- soma_by_gene[!is.na(Mass_Spec_Confirmation)]

sinfo[Mass_Spec_Confirmation != "", Mass_Spec_Confirmation := "SomaLogic"]
sinfo[soma_by_gene, on = .(Aptamer), Mass_Spec_Confirmation := ifelse(Mass_Spec_Confirmation == "", "Emilsson et al. 2018", Mass_Spec_Confirmation)]

# Are the aptamers supported by cis pQTLs?
cis_pQTLs <- fread("analyses/mendelian_randomisation/pqtls/cis_hierarchical_correction.tsv")
cis_pQTLs <- unique(cis_pQTLs[,.(SOMAMER_ID)])
sinfo[, cis_pQTL := ""]
sinfo[cis_pQTLs, on = .(SOMAMER_ID), cis_pQTL := "yes"]

# Order
sinfo <- sinfo[order(-Type)]

#Write out
fwrite(sinfo, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/curated_soma_info.txt")
