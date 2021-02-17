library(data.table)

###############################################
# Get list of PRS associated proteins
###############################################

prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[FDR < 0.05,.(PRS, UniProt, Gene, Beta)])

# Manually add UniProt IDs where there was cross-reactivity
cross_react <- prs_assocs[Gene == "PDE4D"]
cross_react[, c("UniProt", "Gene") := .("P27815", "PDE4A")]
prs_assocs <- rbind(prs_assocs, cross_react)

# Write out
for (this_prs in unique(prs_assocs$PRS)) {
  fwrite(prs_assocs[PRS == this_prs, .(UniProt)], col.names=FALSE, quote=FALSE,
         file=sprintf("analyses/pub/cardiometabolic_proteins/review2/%s_sig_uniprot.txt", this_prs))
}

##############################################################
# Obtain background set of UniProt IDs across all proteins
##############################################################

# Load protein information.
sinfo <- fread("analyses/processed_traits/somalogic_proteins/trait_info.tsv")

# Filter to aptamers passing QC
sinfo <- sinfo[Type == "Protein"]

# Select columns
sinfo <- sinfo[Type == "Protein", .(variable, SOMAMER_ID, Aptamer=SeqId, Target=TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot,
                                    Gene=Gene.Name, Entrez_id=Entrez.Gene.ID, chr, start, end)]

# Fix bad entries (Aptamers for the same target with different/missing gene/uniprot information)
sinfo[Target == "14-3-3 protein family", UniProt := "P61981|Q04917"]
sinfo[Target == "Induced myeloid leukemia cell differentiation protein Mcl-1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("Q07820", "MCL1", "4170", "1", "150547027", "150552214")]
sinfo[Target == "Protein delta homolog 1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("P80370", "DLK1", "8788", "14", "101193202", "101201467")]
sinfo[Target == "Stromal cell-derived factor 1",
  c("UniProt", "Gene", "Entrez_id", "chr", "start", "end") :=
  .("P48061", "CXCL12", "6387", "10", "44865601", "44880545")]

# Obtain list of unique UniProt IDs
background <- data.table(UniProt=unique(unlist(strsplit(sinfo$UniProt, "\\|"))))

# Add entry for PDE4A
background <- rbind(background, data.table(UniProt="P27815"))

# Write out
fwrite(background, col.names=FALSE, quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/all_uniprot.txt")


