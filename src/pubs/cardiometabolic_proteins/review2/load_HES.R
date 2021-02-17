library(data.table)

# Execute all code within environment so we don't override anything in the users session
hes_script_env <- new.env()
with(hes_script_env, {

# Get latest files
hes_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_CaliberEP_.*", full.names=TRUE)
idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
pheno_file <- list.files("data/INTERVAL/HES", pattern="INTERVALdata_[^(p|P)3].*", full.names=TRUE)

# Determine dates
update_date <- data.table(
  file=c(hes_file, idmap_file, pheno_file),
  name=c("HES data", "HES sample identifier linkage file", "HES project phenotype file"),
  date=NA_character_
)
update_date[, date := gsub(".*_", "", file)]
update_date[, date := gsub("\\..*", "", date)]
update_date[, date := ifelse(grepl("[A-Z]", date), 
	as.IDate(date, format="%d%B%Y"),
  as.IDate(date, format="%Y%m%d"))]
update_date[, date := as.IDate(date)] # concatenating two Idates for some reason leads to integer

# Log
update_date[, cat(paste(sprintf("%s last updated: %s\n", name, date), collapse=""))]

# Load HES data and filter to participants with genetic data
hes <- fread(hes_file, colClasses="character")
hes[, identifier := as.integer(identifier)]
idmap <- fread(idmap_file)
hes <- hes[idmap[,.(identifier, IID=Affymetrix_gwasQC_bl)], on=.(identifier)]
hes <- hes[!is.na(IID)]
hes[, identifier := NULL]

# Melt to long
hes <- melt(hes, id.vars="IID", value.name="eventDate")
hes[eventDate == "", eventDate := NA_character_]
hes[, event_code := as.integer(gsub(".*_", "", variable))]
hes[, code_type := gsub(".*_", "", gsub("_d_[0-9]*", "", variable))]
hes[, variable := NULL]
hes[, eventDate := as.IDate(eventDate, format="%d%b%Y")]

# Map event codes to phenotypes
event_codes <- fread("data/INTERVAL/HES/caliber_event_names.csv")
hes <- hes[event_codes, on = .(event_code), nomatch=0]

# get maximum follow-up date
max_followup <- hes[!is.na(eventDate), max(eventDate)]
min_followup <- hes[!is.na(eventDate), min(eventDate)]

# Compute follow-up time
pheno <- fread(pheno_file)
pheno[, attendanceDate := as.IDate(attendanceDate, format="%d%b%Y")]
pheno[, attendanceDate_24m := as.IDate(attendanceDate_24m, format="%d%b%Y")]
pheno[idmap, on = .(identifier), IID := Affymetrix_gwasQC_bl]
pheno <- pheno[!is.na(IID)]
hes[, event := as.integer(!is.na(eventDate))]
hes[is.na(eventDate), eventDate := max_followup]
hes[pheno, on = .(IID), followUp := (eventDate - attendanceDate)/365]
hes[pheno, on = .(IID), followUp_24m := (eventDate - attendanceDate_24m)/365]

# Also age as timescale
hes[pheno, on = .(IID), ageFollowUp := agePulse + followUp]

# Indicate whether event is prevalent or incident
hes[, prevalent := ifelse(followUp < 0, 1L, 0L)]
hes[, prevalent_24m := ifelse(followUp_24m < 0, 1L, 0L)]

# Indicate whether prevalent event was in childhood
hes[, childhood_event := ifelse(ageFollowUp < 18, 1L, 0L)]

# Give code_type more informative label
hes[code_type == "p", code_type := "primary_only"]
hes[code_type == "ps", code_type := "primary+secondary"]

# Extract loaded results out of local script environment
}) # end environment execution

hes <- hes_script_env$hes
rm(hes_script_env)

