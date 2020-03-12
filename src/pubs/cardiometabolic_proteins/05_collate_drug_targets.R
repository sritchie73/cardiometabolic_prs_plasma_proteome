library(data.table)

out_dir <- "analyses/pub/cardiometabolic_proteins"

drugs <- fread(sprintf("%s/drugBank/drug_targets.tsv", out_dir))
trials <- fread(sprintf("%s/drugBank/trial_info.tsv", out_dir))

grs_assocs <- fread(sprintf("%s/all_assocs.tsv", out_dir))
grs_assocs <- grs_assocs[Prot.FDR < 0.05]
grs_assocs <- unique(grs_assocs[,.(PRS, Gene, UniProt)])

# Fix PDE4D entries
grs_assocs <- rbind(grs_assocs, data.table("Chronic Kidney Disease", "PDE4A", "P27815"), use.names=FALSE)
drugs[is.na(PRS), c("PRS", "UniProt", "Gene") := .("Chronic Kidney Disease", "P27815", "PDE4A")]

# Filter to drugs that target the protein, not drugs which are carried, transported,
# or catalyzed by the protein in question
drugs <- drugs[role == "target"]
trials <- trials[DrugBankID %in% drugs$DrugBankID]

# Remove revoked DrugBank entries
revoked <- trials[Status == "revoked", DrugBankID]
drugs <- drugs[!(DrugBankID %in% revoked)]
trials <- trials[!(DrugBankID %in% revoked)]

# Get information for each unique observation
obs <- drugs[, .(PRS, UniProt, Gene, DrugBankID, role, PRSprotPharma)]

# Function to collapse entries in target/enzyme/transporter/carrier list
fixtargets <- function(x) {
  entries <- strsplit(x, "; ")
  sapply(entries, function(entrylist) {
    if (length(entrylist) == 1 && is.na(entrylist)) {
      return(NA)
    }
    tokens <- strsplit(entrylist, " ")
    genes <- sapply(tokens, function(tk) {  tk[length(tk)] })
    types <- sapply(tokens, function(tk) {  if(length(tk) == 1) { "Unknown of" } else { paste(tk[-length(tk)], collapse=" ") } })
    dt <- data.table(type=types, gene=genes)
    dt <- dt[, .(genes=paste(sort(gene), collapse=", ")), by=type]
    dt <- dt[, .(typestring=sort(paste(type, genes)))]
    dt[,paste0(paste(sort(typestring), collapse=". "), ".")]
  })
}

# Collate supp table - each entry is a single drug:
supp <- copy(drugs)
supp[, c("PRS", "UniProt", "Gene", "role") := NULL]
supp <- unique(supp)

# Collapse and organise:
supp[, targets := gsub("^ ", "", fixtargets(targets))]
supp[, enzymes := gsub("^ ", "", fixtargets(enzymes))]
supp[, carriers := gsub("^ ", "", fixtargets(carriers))]
supp[, transporters := gsub("^ ", "", fixtargets(transporters))]

# Collate trial information
phase_order <- data.table(
  Phase = c("Not Available", "0", "1", "1, 2", "2", "2, 3", "3", "4"),
  PhaseNum = 1:8)
trials <- trials[phase_order, on = .(Phase)]

max_phase <- trials[, .SD[which.max(PhaseNum)], by=DrugBankID]
max_phase[, Conditions := gsub("/", " / ", Conditions)]
max_phase <- max_phase[, .(MaxPhase=Phase, 
                           Conditions = paste(Conditions, collapse = " / "),
                           IndicationIDs = paste(IndicationID, collapse=", ")),
          by=.(DrugBankID, Status, Purpose)] # all have unique status and purpose
max_phase <- max_phase[,.(DrugBankID, MaxPhase, Status, Purpose, Indications=Conditions, IndicationIDs)]

supp <- merge(supp, max_phase, by="DrugBankID", all.x=TRUE)

# Collage associated PRS associations
acronym <- function(cvec) {
  words <- strsplit(cvec, " ")
  sapply(words, function(wordlist) {
    first_letters <- sapply(wordlist, function(word) {
      substr(word, 1, 1)
    })
    paste(first_letters, collapse="")
  })
}
assocs <- copy(obs)
assocs[PRSprotPharma == "unknown", PRSprotPharma := "NA"]
assocs[, PRS := acronym(PRS)]
assocs <- assocs[order(Gene)][order(PRS)]
assocs <- assocs[, .(PRS = paste(paste(PRS, collapse=" and "), "PRS")),
              by=.(DrugBankID, UniProt, Gene, PRSprotPharma)]
assocs <- assocs[, .(Genes = paste(Gene, collapse=", "), 
                     UniProts=paste(UniProt, collapse=", "),
                     pharma=paste(PRSprotPharma, collapse=", ")),
                   by=.(DrugBankID, PRS)]
assocs <- assocs[, .(DrugBankID,
  text=sprintf("%s to product of %s [UniProt: %s] [Pharmacological action: %s].", 
               PRS, Genes, UniProts, pharma))]
assocs <- assocs[, .(PRS_association = paste(sort(text), collapse=" ")), by=DrugBankID]

supp <- merge(assocs, supp, by="DrugBankID", all.y=TRUE)

# Curate list of indications relevant to each PRS
conds <- trials[obs[, .(DrugBankID, UniProt, Gene, PRS)], on=.(DrugBankID), allow.cartesian=TRUE, nomatch=0]
conds <- conds[, .(Indication=strsplit(Conditions, "/")[[1]]), by=.(PRS, UniProt, Gene, DrugBankID, IndicationID)]

t2d_conds <- conds[PRS == "Type 2 Diabetes"]
ckd_conds <- conds[PRS == "Chronic Kidney Disease"]
cad_conds <- conds[PRS == "Coronary Artery Disease"]

filt <- t2d_conds[grepl("diabetes", Indication, ignore.case=TRUE)]
filt <- filt[!grepl("type 1", Indication, ignore.case=TRUE)]
filt <- filt[!grepl("Type1", Indication, ignore.case=TRUE)]
filt <- filt[!grepl("insipidus", Indication, ignore.case=TRUE)]
filt <- filt[!grepl("gestation", Indication, ignore.case=TRUE)]
filt <- filt[!grepl("pregnan", Indication, ignore.case=TRUE)]
filt <- filt[!grepl("autoimmune", Indication, ignore.case=TRUE)]

filt <- rbind(filt, t2d_conds[grepl("insulin", Indication, ignore.case=TRUE)])
filt <- filt[!grepl("insulinoma", Indication, ignore.case=TRUE)]
filt <- filt[!grepl("Type 1", Indication, ignore.case=TRUE)]
filt <- filt[!grepl("Insulin-Like", Indication, ignore.case=TRUE)]

filt <- rbind(filt, t2d_conds[grepl("glucose", Indication, ignore.case=TRUE)])
filt <- rbind(filt, t2d_conds[grepl("glyca?emia", Indication, ignore.case=TRUE)])
filt <- rbind(filt, t2d_conds[grepl("insulina?emia", Indication, ignore.case=TRUE)])

filt <- rbind(filt, ckd_conds[grepl("kidney", Indication, ignore.case=TRUE)])
filt <- rbind(filt, ckd_conds[grepl("renal", Indication, ignore.case=TRUE)])
filt <- rbind(filt, ckd_conds[grepl("neph", Indication, ignore.case=TRUE)])

filt <- rbind(filt, cad_conds[grepl("cardio", Indication, ignore.case=TRUE)])
filt <- filt[!grepl("Cardiotoxicity", Indication)]
filt <- rbind(filt, cad_conds[grepl("coronary", Indication, ignore.case=TRUE)])
filt <- rbind(filt, cad_conds[grepl("vascular", Indication, ignore.case=TRUE)])
filt <- filt[!grepl("Accident", Indication)]
filt <- filt[!grepl("Macul", Indication)]
filt <- rbind(filt, cad_conds[grepl("heart", Indication, ignore.case=TRUE)])
filt <- filt[!grepl("Congenital", Indication)]
filt <- rbind(filt, cad_conds[grepl("arter", Indication, ignore.case=TRUE)])
filt <- rbind(filt, cad_conds[grepl("vessel", Indication, ignore.case=TRUE)])

filt <- unique(filt)
filt <- filt[,Indication := NULL]

# Get max phase per condition per drug
prs_trials <- trials[filt, on = .(DrugBankID, IndicationID)]
prs_trials <- prs_trials[, .SD[which.max(PhaseNum)], by=.(PRS, UniProt, Gene, DrugBankID, IndicationID)]
prs_trials <- prs_trials[order(-PhaseNum)]
prs_trials[, PRS_trials := sprintf("[Max phase: %s, Purpose: %s, Status: %s, Indication: %s] %s.",
                                   Phase, Purpose, Status, IndicationID, gsub("/", " / ", Conditions))]
prs_trials <- prs_trials[, .(PRS_trials = paste(PRS_trials, collapse=" ")), by=.(DrugBankID)]

supp <- merge(supp, prs_trials, by="DrugBankID", all.x=TRUE)

# Order and select
group_lvls <- data.table(
  groups = c("experimental", "experimental; investigational", "investigational",
             "approved; withdrawn", "approved; investigational; withdrawn",
             "approved; experimental", "approved; experimental; investigational",
             "approved; investigational", "approved; investigational; nutraceutical",
             "approved; nutraceutical", "approved"),
  gnum = 1:11)

supp[group_lvls, on = .(groups), group_num := gnum]
supp[, groups := gsub(";", ",", groups)]
supp[phase_order, on = .(MaxPhase=Phase), PhaseNum := PhaseNum]
supp <- supp[order(!is.na(PRS_trials))][order(-group_num)][order(-PhaseNum)]

supp <- supp[,.(DrugBankID, PRS_association, name, description, type, groups, indication,
                pharmacodynamics, mechanism, pharmaAction, targets, enzymes, transporters,
                carriers, MaxPhase, Status, Purpose, Indications, IndicationIDs, PRS_trials)]

# Write out 
fwrite(supp, sep="\t", quote=FALSE, file=sprintf("%s/all_drug_targets.tsv", out_dir))

# Summarise drugs targetting each protein
main_prot <- obs[, .(N_drugs=.N, N_pharma=sum(PRSprotPharma == "yes")), by= .(PRS, UniProt, Gene)]
approved <- supp[groups %like% "approved" & !(groups %like% "withdrawn"), .(DrugBankID)]
approved <- obs[approved, on = .(DrugBankID)]
approved <- approved[, .(N_approved=.N, N_approved_pharma=sum(PRSprotPharma == "yes")), by= .(PRS, UniProt, Gene)]
main_prot <- merge(main_prot, approved, by=c("PRS", "UniProt", "Gene"), all.x=TRUE)
main_prot[is.na(N_approved), c("N_approved", "N_approved_pharma") := 0]
main_prot[filt, on = .(PRS, UniProt, Gene), relevant_trial := TRUE]
main_prot[!filt, on = .(PRS, UniProt, Gene), relevant_trial := FALSE]
main_prot <- main_prot[order(Gene)][order(PRS)]

fwrite(main_prot, sep="\t", quote=FALSE, file=sprintf("%s/prots_drug_targets.tsv", out_dir))

# Get key list of drugs through which their pharmacological action is through a PRS associated target
main_drugs <- obs[PRSprotPharma == "yes"]
main_drugs <- supp[main_drugs, on =.(DrugBankID)]
main_drugs <- main_drugs[, .(PRS, UniProt, Gene, DrugBankID, name, type, 
                             groups, indication, mechanism, MaxPhase, PRS_trials)]

fwrite(main_drugs, sep="\t", quote=FALSE, file=sprintf("%s/main_drugs.tsv", out_dir))

# Get full list of drugs to curate a list of drugs through which the mechanism
# is opposed to the PRS


