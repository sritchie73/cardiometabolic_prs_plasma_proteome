library(XML)
library(data.table)
library(foreach)

out_dir <- "analyses/pub/cardiometabolic_proteins/drugBank/"
dir.create(out_dir, showWarnings=FALSE)

# This takes a while to load.
# See https://docs.drugbankplus.com/xml for tree structure
data <- xmlTreeParse("data/DrugBank/DrugBank_full_db_20190624.xml")
root <- xmlRoot(data)

# Uniprot IDs of proteins to search for
grs_assocs <- fread("analyses/pub/cardiometabolic_proteins/all_assocs.tsv")
grs_assocs <- grs_assocs[Prot.FDR < 0.05]
prots <- unique(unlist(strsplit(grs_assocs$UniProt, "\\|")))
prots <- c(prots, "P27815") # add PDE4A

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
      actions <- ifelse(length(actions) == 0, NA_character_, paste(actions, collapse=", "))

      tryCatch({
        UniProt <- xmlAttrs(targetNode[["polypeptide"]])[["id"]]
      }, error=function(e) {
        UniProt <- NA_character_
      })
      gene <- xmlValue(targetNode[["polypeptide"]][["gene-name"]])
      if (length(gene) == 0) gene <- NA_character_

      data.table(DrugBankID = id, UniProt, gene, type=targetType, organism, pharmaAction, actions)
    }
  }
}

# Combine information:
combined <- grs_assocs[,.(UniProt=strsplit(UniProt, ",")[[1]]), by=.(PRS, Gene)]
combined <- combined[drugTargets, on = .(UniProt=target), nomatch=0, .(PRS, UniProt, Gene, DrugBankID)] 
combined <- combined[drugInfo, on = .(DrugBankID)]

# What role does the PRS-associated protein play? Is it a target, enzyme, carrier, or transporter?
# Does the drug target that protein in humans or some other organism?
combined[targetInfo, on = .(DrugBankID, UniProt), role := gsub("s$", "", i.type)]
combined[targetInfo, on = .(DrugBankID, UniProt), role_organism := organism] # check PRS-associated protein is human target
combined[is.na(role_organism) | role_organism != "Humans", role := paste0(role, " (", role_organism, ")")]
combined[, role_organism := NULL]

# What are the pharmacological targets of the drug?
capitalize <- function(x) { paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))) }
pharmaInfo <- targetInfo[pharmaAction == "yes", .(DrugBankID, organism, actions, gene)]
pharmaInfo[, text := gene]
pharmaInfo[!is.na(actions), text := paste(capitalize(actions), "of", gene)]
pharmaInfo[organism != "Humans", text := paste(text, "in", organism)]
pharmaInfo <- pharmaInfo[, .(text=paste(text, collapse="; ")), by=.(DrugBankID)]
combined[pharmaInfo, on = .(DrugBankID), pharmaAction := text]

# Add other target information:
targetText <- copy(targetInfo)
targetText[, text := gene]
targetText[!is.na(actions), text := paste(capitalize(actions), "of", gene)]
targetText[organism != "Humans", text := paste(text, "in", organism)]
targetText <- targetText[pharmaAction != "yes", .(text=paste(text, collapse="; ")), by=.(DrugBankID, type)]
combined[targetText[type == "targets"], on = .(DrugBankID), targets := text]
combined[targetText[type == "enzymes"], on = .(DrugBankID), enzymes := text]
combined[targetText[type == "transporters"], on = .(DrugBankID), transporters := text]
combined[targetText[type == "carriers"], on = .(DrugBankID), carriers := text]

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

# Output
fwrite(combined, sep="\t", quote=FALSE, na="NA", file=sprintf("%s/drug_targets.tsv", out_dir))

# Look up clinical trial information from each drugs corresponding page
# (data only provided in XML file to paying customers)
trialInfo <- foreach(drugID = unique(na.omit(combined$DrugBankID)), .combine=rbind) %do% {
  system(sprintf("wget -O %s/%s.html https://www.drugbank.ca/drugs/%s", out_dir, drugID, drugID))
  html <- htmlTreeParse(sprintf("%s/%s.html", out_dir, drugID))
  system(sprintf("rm %s/%s.html", out_dir, drugID))

  # check for revocation status
  if (grepl("revoked", xmlValue(xmlRoot(html)[["body"]][["main"]][[1]][[1]][[1]]))) {
    dt <- data.table(Phase = NA, status = "revoked", Purpose = NA, 
                     Conditions = NA, Count = NA, IndicationID = NA,
                     DrugBankID = drugID)
    return(dt)
  }

  # Otherwise extract main content on webpage
  content <- xmlRoot(html)[["body"]][["main"]][["div"]][[4]]

  # Find clinical trial heading, if available
  tIdx <- NULL
  for (nodeIdx in seq_along(content)) {
    if (names(content)[[nodeIdx]] == "h2") {
      if (xmlAttrs(content[[nodeIdx]])[["id"]] == "clinical-trials") {
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
  if(!is.data.table(tbl)) { stop(drugID) }
  setnames(tbl, c("Phase", "Status", "Purpose", "Conditions", "Count", "IndicationID"))

  tbl[, DrugBankID := drugID]
  tbl
}

fwrite(trialInfo, sep="\t", quote=FALSE, file=sprintf("%s/trial_info.tsv", out_dir))
