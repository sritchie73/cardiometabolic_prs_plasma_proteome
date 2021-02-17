library(XML)
library(data.table)
library(foreach)
library(doMC)

# Set up parallel environment
if (!exists("ncores")) {
  ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
  ncores <- as.integer(ncores)
  if(is.na(ncores)) ncores <- 1
}
registerDoMC(ncores)
setDTthreads(ncores)

# Set up output directory
out_dir <- "analyses/pub/cardiometabolic_proteins/review2/drugbank"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# This takes a while to load.
# See https://docs.drugbankplus.com/xml for tree structure
data <- xmlTreeParse("data/DrugBank/DrugBank_full_db_v5.17_20200702.xml")
root <- xmlRoot(data)

# Load FDR significant associations
prs_assocs <- fread("analyses/pub/cardiometabolic_proteins/review2/all_prs_to_prot_assocs.txt")
prs_assocs <- unique(prs_assocs[,.(PRS, Target, UniProt, Gene, PRS.Beta=Beta, PRS.L95=L95, PRS.U95=U95, PRS.P=P, PRS.FDR=FDR)])
prs_assocs <- prs_assocs[PRS.FDR < 0.05]

# Add in PDE4A
pde4a <- prs_assocs[Gene == "PDE4D"]
pde4a[, c("Target", "UniProt", "Gene") := .("cAMP-specific 3',5'-cyclic phosphodiesterase 4A", "P27815", "PDE4A")]
prs_assocs <- rbind(prs_assocs, pde4a)

# Uniprot IDs of proteins to search for
prots <- unique(unlist(strsplit(prs_assocs$UniProt, "\\|")))

# Iterate through tree finding drugs that target the list of proteins
drugTargets <- foreach(drugIdx = seq_along(root), .combine=rbind) %do% {
  drugNode <- root[[drugIdx]]

  # First get a vector of UniProt ids for the drug's targets, enzymes, carriers, and transporters:
  targets <- foreach(targetType = c("targets", "transporters", "enzymes", "carriers"), .combine=rbind) %do% {
    foreach(targetIdx = seq_along(drugNode[[targetType]]), .combine=rbind, .errorhandling="remove") %do% {
      # Sometimes fails where the polypeptide node is NULL/does not exist.
      xmlAttrs(drugNode[[targetType]][[targetIdx]][["polypeptide"]])[["id"]]
    }
  }

  # Find which target(s) of interest correspond to the drug
  if (any(prots %in% targets)) {
    data.table(nodeIdx = drugIdx, 
               DrugBankID = xmlValue(drugNode[["drugbank-id"]]), 
               target = intersect(prots, targets))
  }
}

# Get basic information about each drug
drugInfo <- foreach(drugIdx = unique(drugTargets$nodeIdx), .combine=rbind) %do% {
  drugNode <- root[[drugIdx]]
  id <- xmlValue(drugNode[["drugbank-id"]])
  name <- xmlValue(drugNode[["name"]])
  type <- xmlAttrs(drugNode)[["type"]]

  xmlValueMaybe <- function(...) {
    value <- XML::xmlValue(...)
    if (length(value) == 0) {
      return(NA_character_)
    } else {
      return(gsub("\r\n\r\n", " ", value))
    }
  }

  description <- xmlValueMaybe(drugNode[["description"]])
  indication <- xmlValueMaybe(drugNode[["indication"]])
  pharmacodynamics <- xmlValueMaybe(drugNode[["pharmacodynamics"]])
  mechanism <- xmlValueMaybe(drugNode[["mechanism-of-action"]])

  groups <- foreach(gIdx = seq_along(drugNode[["groups"]]), .combine=c) %do% {
    xmlValue(drugNode[["groups"]][[gIdx]])
  }
  groups <- groups[groups != "vet_approved"] # don't care about animal use
  groups <- paste(groups, collapse="; ")

  organisms <- foreach(oIdx = seq_along(drugNode[["affected-organisms"]]), .combine=c) %do% {
    xmlValue(drugNode[["affected-organisms"]][[oIdx]])
  }
  organisms <- ifelse(length(organisms) > 0, paste(organisms, collapse="; "), NA_character_)

  biotech <- foreach(bIdx = seq_along(drugNode[["biotech-categories"]]), .combine=c) %do% {
    xmlValue(drugNode[["biotech-categories"]][[bIdx]])
  }
  biotech <- ifelse(length(biotech) > 0, paste(biotech, collapse="; "), NA_character_)

  ChEMBLID <- foreach(eIdx = seq_along(drugNode[["external-identifiers"]]), .combine=c) %do% {
    node <- drugNode[["external-identifiers"]][[eIdx]]
    if (xmlValue(node[["resource"]]) == "ChEMBL") {
      return(xmlValue(node[["identifier"]]))
    }
  }
  if (is.null(ChEMBLID)) ChEMBLID <- NA_character_

  data.table(DrugBankID = id, name, description, type, groups, organisms, 
             biotech, indication, pharmacodynamics, mechanism, ChEMBLID)
}

# Functions to aid conversion of data to text
capitalize <- function(x) { 
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))) 
}

wordlist <- function(x) {
  last_two <- x[(length(x) - 1):length(x)]
  rest <- x[-((length(x) - 1):length(x))]
  if (length(rest) > 0) {
		paste0(paste(rest, collapse=", "), ", ", paste(last_two, collapse=", and "))
  } else {
    paste(last_two, collapse=" and ")
  }
}

pluralfy <- function(x, n) {
  if (n == 1) {
   gsub("s$", "", x)
  } else {
   x
  }
}

# Get information about all targets for each drug
targetInfo <- foreach(drugIdx = unique(drugTargets$nodeIdx), .combine=rbind) %do% {
  drugNode <- root[[drugIdx]]
  id <- xmlValue(drugNode[["drugbank-id"]])

  foreach(targetType = c("targets", "transporters", "enzymes", "carriers"), .combine=rbind) %do% {
    foreach(targetIdx = seq_along(drugNode[[targetType]]), .combine=rbind) %do% {
      targetNode <- drugNode[[targetType]][[targetIdx]]
      organism <- xmlValue(targetNode[["organism"]])
      if (length(organism) == 0) organism <- NA_character_
      pharmaAction <- xmlValue(targetNode[["known-action"]]) # is the pharmacological action of this drug due to this target?
      actions <- foreach(aIdx = seq_along(targetNode[["actions"]]), .combine=c) %do% {
        xmlValue(targetNode[["actions"]][[aIdx]])
      }
      actions <- ifelse(length(actions) == 0, "unclassified", wordlist(actions))

      tryCatch({
        UniProt <- xmlAttrs(targetNode[["polypeptide"]])[["id"]]
      }, error=function(e) {
        UniProt <- NA_character_
      })
      gene <- xmlValue(targetNode[["polypeptide"]][["gene-name"]])

      if (length(gene) == 0 || is.na(gene)) {
        gene <- if (UniProt == "P16220") { "CREB1" 
        } else if (UniProt == "P19099") { "CYP11B2" 
        } else if (UniProt == "Q06432") { "CACNG1" 
        } else if (UniProt == "P01619") { "IGKV3-20" 
        } else if (UniProt == "Q9UPY5") { "SLC7A11" 
        } else if (UniProt == "Q07343") { "PDE4B" 
        } else if (UniProt == "P27815") { "PDE4A" 
        } else { NA_character_ }
      }

      data.table(DrugBankID = id, UniProt, gene, type=targetType, organism, pharmaAction, actions)
    }
  }
}

# Filter drug targets to those in humans
targetInfo <- targetInfo[organism == "Humans"]

# Collate the few multiple entries
targetInfo <- unique(targetInfo)[, .(actions=wordlist(actions)), by = .(DrugBankID, UniProt, gene, type, pharmaAction)]

# Split into wide format by drug
targetWide <- unique(targetInfo[,.(DrugBankID)])

pharmaAct <- targetInfo[pharmaAction == "yes", .(text=paste(capitalize(actions), "of", pluralfy(type, .N), wordlist(sort(gene)))), by = .(DrugBankID, type, actions)]
pharmaAct <- pharmaAct[order(actions)][order(-type)][, .(text=paste0(paste(text, collapse="; "), ".")), by=.(DrugBankID)]
targetWide[pharmaAct, on = .(DrugBankID), pharmaAct := i.text]

otherTargets <- targetInfo[pharmaAction != "yes" & type == "targets", .(text=paste0(capitalize(actions), " of ", wordlist(sort(gene)))), by = .(DrugBankID, actions)]
otherTargets <- otherTargets[order(actions)][, .(text=paste0(paste(text, collapse="; "), ".")), by=.(DrugBankID)]
targetWide[otherTargets, on = .(DrugBankID), Targets := i.text]

otherEnzymes <- targetInfo[pharmaAction != "yes" & type == "enzymes", .(text=paste0(capitalize(actions), " of ", wordlist(sort(gene)))), by = .(DrugBankID, actions)]
otherEnzymes <- otherEnzymes[order(actions)][, .(text=paste0(paste(text, collapse="; "), ".")), by=.(DrugBankID)]
targetWide[otherEnzymes, on = .(DrugBankID), Enzymes := i.text]

otherCarriers <- targetInfo[pharmaAction != "yes" & type == "carriers", .(text=paste0(capitalize(actions), " of ", wordlist(sort(gene)))), by = .(DrugBankID, actions)]
otherCarriers <- otherCarriers[order(actions)][, .(text=paste0(paste(text, collapse="; "), ".")), by=.(DrugBankID)]
targetWide[otherCarriers, on = .(DrugBankID), Carriers := i.text]

otherTransporters <- targetInfo[pharmaAction != "yes" & type == "transporters", .(text=paste0(capitalize(actions), " of ", wordlist(sort(gene)))), by = .(DrugBankID, actions)]
otherTransporters <- otherTransporters[order(actions)][, .(text=paste0(paste(text, collapse="; "), ".")), by=.(DrugBankID)]
targetWide[otherTransporters, on = .(DrugBankID), Transporters := i.text]

# Combine information for supp table
combined <- prs_assocs[,.(UniProt=strsplit(UniProt, ",")[[1]]), by=.(PRS, Target, UniProt, Gene)]
combined <- combined[drugTargets, on = .(UniProt=target), nomatch=0, .(PRS, UniProt, Gene, DrugBankID), allow.cartesian=TRUE] 
combined <- combined[drugInfo, on = .(DrugBankID)]
combined <- combined[targetWide, on = .(DrugBankID)]

# What's the pharamcologic role of the PRS-associated protein?
combined[targetInfo, on = .(DrugBankID, UniProt), PRSprotPharma := i.pharmaAction]

# Strip newline characters from descriptive fields
combined[, description := gsub("\n", "    ", description)]
combined[, indication := gsub("\n", "    ", indication)]
combined[, pharmacodynamics := gsub("\n", "    ", pharmacodynamics)]
combined[, mechanism := gsub("\n", "    ", mechanism)]
combined[, description := gsub("\r", "", description)]
combined[, indication := gsub("\r", "", indication)]
combined[, pharmacodynamics := gsub("\r", "", pharmacodynamics)]
combined[, mechanism := gsub("\r", "", mechanism)]

# Look up clinical trial information from each drugs corresponding page
# (data only provided in XML file to paying customers)
trialInfo <- foreach(drugID = unique(na.omit(combined$DrugBankID)), .combine=rbind) %do% {
  system(sprintf("wget -O %s/%s.html https://www.drugbank.ca/drugs/%s", out_dir, drugID, drugID))
  html <- htmlTreeParse(sprintf("%s/%s.html", out_dir, drugID))
  system(sprintf("rm %s/%s.html", out_dir, drugID))

  # check for revocation status
  if (grepl("revoked", xmlValue(xmlRoot(html)[["body"]][["main"]][[1]][[1]][[1]]))) {
    dt <- data.table(Phase = NA, Status = "revoked", Purpose = NA,
                     Conditions = NA, Count = NA, IndicationID = NA,
                     DrugBankID = drugID)
    return(dt)
  }

  # Otherwise extract main content on webpage
  content <- xmlRoot(html)[["body"]][["main"]][["div"]][["div"]][[2]][[2]]

  # Find clinical trial heading, if available
  tIdx <- NULL
  for (nodeIdx in seq_along(content)) {
    if (names(content)[[nodeIdx]] == "h2") {
      if (xmlAttrs(content[[nodeIdx]])[["id"]] == "clinical-trials-header") {
         tIdx <- nodeIdx + 1
      }
    }
  }

  # Extract table
  tblNode <- content[[tIdx]][["dd"]][["table"]]
  if (is.null(tblNode)) {
    return(NULL)
  }

  tbl <- foreach(rIdx = seq_along(tblNode[["tbody"]]), .combine=rbind) %do% {
    l <- foreach(cIdx = seq_along(tblNode[["tbody"]][[rIdx]])) %do% {
      xmlValue(tblNode[["tbody"]][[rIdx]][[cIdx]])
    }
    dt <- as.data.table(l)

    id <- xmlAttrs(tblNode[["tbody"]][[rIdx]][[4]][[1]])[["href"]]
    id <- gsub("/indications/", "", id)

    dt[, id := id]

    return(dt)
  }
  if (is.null(tbl)) {
    return(NULL)
  }

  setnames(tbl, c("Phase", "Status", "Purpose", "Conditions", "Count", "IndicationID"))
  tbl[, DrugBankID := drugID]
  tbl
}

# Drop revoked drug(s)
combined <- combined[!(DrugBankID %in% trialInfo[Status == "revoked", DrugBankID])]
trialInfo <- trialInfo[Status != "revoked"]

# Get approved indications or maximum phase clinic trial indications
phase_order <- data.table(
  Phase = c("Not Available", "0", "1", "1, 2", "2", "2, 3", "3", "4"),
  PhaseNum = 1:8)
trialInfo <- trialInfo[phase_order, on = .(Phase)]

max_phase <- trialInfo[, .SD[PhaseNum == max(PhaseNum)], by=DrugBankID]
max_phase[, Conditions := gsub("/", " / ", Conditions)]
max_phase <- max_phase[, .(MaxPhase=unique(Phase), 
                           Status = paste(Status, collapse = " | "),
                           Purpose = paste(Purpose, collapse = " | "),
                           Indications = paste(Conditions, collapse = " | "),
                           IndicationIDs = paste(IndicationID, collapse=" | ")),
          by=.(DrugBankID)]

# Add to supp table
combined <- merge(combined, max_phase, by = "DrugBankID", all.x = TRUE)

# Curate list of indications relevant to each PRS
conds <- trialInfo[combined[, .(DrugBankID, UniProt, Gene, PRS)], on=.(DrugBankID), allow.cartesian=TRUE, nomatch=0]
conds <- conds[, .(Indication=strsplit(Conditions, "/")[[1]]), by=.(PRS, UniProt, Gene, DrugBankID, IndicationID)]

t2d_conds <- conds[PRS == "T2D_PRS"]
ckd_conds <- conds[PRS == "CKD_PRS"]
cad_conds <- conds[PRS == "CAD_PRS"]
is_conds <- conds[PRS == "IS_PRS"]

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

filt <- rbind(filt, is_conds[grepl("ischa?emic", Indication, ignore.case=TRUE)])

filt <- unique(filt)
filt <- filt[,Indication := NULL]

# Get max phase per condition per drug
prs_trials <- trialInfo[filt, on = .(DrugBankID, IndicationID)]
prs_trials <- prs_trials[, .SD[PhaseNum == max(PhaseNum)], by=.(PRS, UniProt, Gene, DrugBankID, IndicationID)]
prs_trials <- prs_trials[order(-PhaseNum)]
prs_trials[, PRS_trials := sprintf("[Max phase: %s, Purpose: %s, Status: %s, Indication: %s] %s.",
                                   Phase, Purpose, Status, IndicationID, gsub("/", " / ", Conditions))]
prs_trials <- prs_trials[, .(PRS_trials = paste(PRS_trials, collapse=" ")), by=.(PRS, DrugBankID)]

# Add to supp table
combined <- merge(combined, prs_trials, by = c("PRS","DrugBankID"), all.x = TRUE)

# Get list of drugs whose effect on the protein is consistent with the PRS association
prs_target_effects <- prs_assocs[,.(PRS, UniProt, PRS.Beta)]
prs_target_effects <- prs_target_effects[targetInfo, on = .(UniProt), nomatch=0]

prs_target_effects[, drug_effect_direction := "unknown"]
prs_target_effects[actions %in% c("agonist", "binder and potentiator", "inducer", "ligand", "product of"), drug_effect_direction := "positive"]
prs_target_effects[actions %in% c("antagonist", "antagonist and inhibitory allosteric modulator", "inhibitor",
                                  "inhibitor and binder", "substrate"), drug_effect_direction := "negative"]
prs_target_effects[actions %in% c("binder", "carrier"), drug_effect_direction := "effected by"]
prs_target_effects[, counterPRS := FALSE]
prs_target_effects[PRS.Beta < 0 & drug_effect_direction == "positive", counterPRS := TRUE]
prs_target_effects[PRS.Beta > 0 & drug_effect_direction == "negative", counterPRS := TRUE]

# Collate PRS target relevant information
prsTargetInfo <- prs_assocs[,.(PRS, Protein=Target, UniProt, Gene)]
prsTargetInfo <- prsTargetInfo[combined[,.(PRS, Gene, DrugBankID, DrugName=name)], on = .(PRS, Gene)]
prsTargetInfo <- prsTargetInfo[targetInfo, on = .(DrugBankID, UniProt), nomatch=0]
prsTargetInfo[, TargetText := sprintf("%s of %s %s.", capitalize(actions), pluralfy(type, 1), Gene)]
prsTargetInfo[prs_target_effects, on = .(PRS, UniProt, DrugBankID), c("effectDirection", "drugOpposePRS") := .(drug_effect_direction, counterPRS)]

# Fix duplicates
prsTargetInfo <- prsTargetInfo[, .(TargetText=paste(TargetText, collapse=" "),
                                   effectDirection=paste(unique(effectDirection), collapse="| "),
                                   pharmaAction = paste(unique(pharmaAction), collapse="| "),
                                   PRS = unique(PRS),
                                   drugOpposePRS=all(drugOpposePRS)),
                                by=.(PRS, Protein, UniProt, Gene, DrugBankID, DrugName)]
prsTargetInfo[DrugName == "NADH" & Gene == "ADH4", TargetText := "Substrate of enzyme ADH4."]

# Fix other problems
prsTargetInfo[grepl("Product of of", TargetText), TargetText := gsub("of of", "of", TargetText)]
prsTargetInfo[grepl("Carrier of carrier", TargetText), TargetText := gsub("Carrier of carrier", "Carried by", TargetText)]

# Add extended drug information
prsTargetInfo <- prsTargetInfo[combined, on = .(PRS, UniProt, Gene, DrugBankID, DrugName=name)]

# Write out supp table
fwrite(prsTargetInfo, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/PRS_to_drug_targets.txt")

# Collate drug target stats
prot_stats <- unique(prsTargetInfo[,.(Protein, UniProt, PRS)])
n_drugs <- prsTargetInfo[,.(N=length(unique(DrugName))), by=UniProt]
prot_stats[n_drugs, on=.(UniProt), N_drugs := i.N]
n_appr_drugs <- prsTargetInfo[groups %like% "approved" & !(groups %like% "withdrawn"), .(N=length(unique(DrugName))), by=UniProt]
prot_stats[n_appr_drugs, on=.(UniProt), N_approved := i.N]
prot_stats[, prs_rt := "No"]
prot_stats[prsTargetInfo[!is.na(PRS_trials)], on = .(UniProt, PRS), prs_rt := "Yes"]

fwrite(prot_stats, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/drug_target_summary.txt")






