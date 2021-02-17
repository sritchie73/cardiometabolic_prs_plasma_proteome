library(data.table)

source("src/utilities/load_HES.R")
source("src/pubs/cardiometabolic_proteins/review2/cardiometabolic_events.R")

# Get list of events
prev_tab <- fread("data/INTERVAL/HES/caliber_event_names.csv")

# Get lists of ICD10 codes
gh <- list.files(path="data/CALIBER_github/secondary_care", pattern="ICD_*", full.names=TRUE)
icd10 <- lapply(gh, fread)
names(icd10) <- gh
icd10 <- rbindlist(idcol="path", icd10)
icd10[, short_name := gsub("ICD_", "", gsub(".csv", "", basename(path)))]
icd10 <- icd10[, .(ICD10 = paste(sort(ICD10code), collapse=", ")), by=short_name]
icd10[short_name == "2ry_polycythaemia", short_name := "sec_polycythaemia"]

# Curate list of events
prev_tab <- prev_tab[icd10, on = .(phenotype_short=short_name)]
prev_tab[, c("phenotype_short", "event_code") := NULL]
prev_tab[, cardiometabolic := phenotype %in% cardiometabolic]

# Get all prevalent events
prev <- unique(hes[prevalent == 1, .(IID, phenotype)])

idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
idmap <- fread(idmap_file)

prev <- prev[IID %in% idmap$Affymetrix_gwasQC_bl]
prev[idmap, on = .(IID=Affymetrix_gwasQC_bl), has_somalogic := !is.na(soma4000_gwasQC_bl)]

# Tabulate numbers in whole of INTERVAL and in subset with somalogic data
prev_n <- prev[,.(N_all=.N, N_soma=sum(has_somalogic)), by=.(phenotype)]
prev_tab <- merge(prev_tab, prev_n, by="phenotype", all.x=TRUE)
prev_tab[is.na(N_all), N_all := 0]
prev_tab[is.na(N_soma), N_soma := 0]

# Get totals
totals <- idmap[!is.na(Affymetrix_gwasQC_bl), .(Total_all = .N, Total_soma = sum(!is.na(soma4000_gwasQC_bl)))]
total_excl <- data.table(
  Excl_all = prev[phenotype %in% cardiometabolic, length(unique(IID))],
  Excl_soma = prev[phenotype %in% cardiometabolic & has_somalogic, length(unique(IID))]
)
total_kept <- data.table(
  Kept_all = totals$Total_all - total_excl$Excl_all,
  Kept_soma = totals$Total_soma - total_excl$Excl_soma
)
prev_any <- data.table(
  Prev_all = prev[, length(unique(IID))],
  Prev_soma = prev[(has_somalogic), length(unique(IID))]
)

# Sort events by number and category
prev_tab <- prev_tab[order(-N_all)][order(-cardiometabolic)]

# Split by category
prev_cardio <- prev_tab[(cardiometabolic)]
prev_other <- prev_tab[!(cardiometabolic)]

# Add percentages for events excluded
prev_cardio[, N_all := sprintf("%s (%s%%)", format(N_all, big.mark=","), round(N_all/total_excl$Excl_all, digits=3)*100)]
prev_cardio[, N_soma := sprintf("%s (%s%%)", format(N_soma, big.mark=","), round(N_soma/total_excl$Excl_soma, digits=3)*100)]

prev_other[, N_all := sprintf("%s (%s%%)", format(N_all, big.mark=","), round(N_all/prev_any$Prev_all, digits=3)*100)]
prev_other[, N_soma := sprintf("%s (%s%%)", format(N_soma, big.mark=","), round(N_soma/prev_any$Prev_soma, digits=3)*100)]

# Curate totals
total_cardio <- data.table(
  phenotype = c("Total participants:", "Prevalent cardiometabolic disease:", "No cardiometabolic disease:"),
  N_all = c(format(totals$Total_all, big.mark=","), 
            sprintf("%s (%s%%)", format(total_excl$Excl_all, big.mark=","), round(total_excl$Excl_all/totals$Total_all, digits=3)*100), 
            sprintf("%s (%s%%)", format(total_kept$Kept_all, big.mark=","), round(total_kept$Kept_all/totals$Total_all, digits=3)*100)),
  N_soma = c(format(totals$Total_soma, big.mark=","), 
            sprintf("%s (%s%%)", format(total_excl$Excl_soma, big.mark=","), round(total_excl$Excl_soma/totals$Total_soma, digits=3)*100), 
            sprintf("%s (%s%%)", format(total_kept$Kept_soma, big.mark=","), round(total_kept$Kept_soma/totals$Total_soma, digits=3)*100))
)

total_other <- data.table(
  phenotype = "Any prevalent disease:",
  N_all = sprintf("%s (%s%%)", format(prev_any$Prev_all, big.mark=","), round(prev_any$Prev_all/totals$Total_all, digits=3)*100),
  N_soma = sprintf("%s (%s%%)", format(prev_any$Prev_soma, big.mark=","), round(prev_any$Prev_soma/totals$Total_soma, digits=3)*100)
)
  
# Build output table
prev_tab <- rbind(fill=TRUE, total_cardio, prev_cardio, total_other, prev_other)
prev_tab <- prev_tab[,.(phenotype, ICD10, cardiometabolic, N_all, N_soma)]

fwrite(prev_tab, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review2/prevalent_tabulated.txt")




