library(data.table)
library(lubridate)

# Execute all code within environment so we don't override anything in the users session
hes_script_env <- new.env()
with(hes_script_env, {

# Get latest files
hes_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_CaliberEP_.*", full.names=TRUE)
covid_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_Covid_.*", full.names=TRUE)
idmap_file <- list.files("data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[^(p|P)3].*", full.names=TRUE)
pheno_file <- list.files("data/INTERVAL/HES", pattern="INTERVALdata_[^(p|P)3].*", full.names=TRUE)

# Determine dates
update_date <- data.table(
  file=c(hes_file, covid_file, idmap_file, pheno_file),
  name=c("HES data", "COVID test data", "HES sample identifier linkage file", "HES project phenotype file"),
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
hes[pheno, on = .(IID), followUp := time_length(as.Date(eventDate) - as.Date(attendanceDate), unit="year")]
hes[pheno, on = .(IID), followUp_24m := time_length(as.Date(eventDate) - as.Date(attendanceDate_24m), unit="year")]

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

# Load covid19 date and determine whether the incident event occurs post-covid19 positive test
covid <- fread(covid_file)
covid[idmap, on = .(identifier), IID := Affymetrix_gwasQC_bl]
covid <- melt(covid, id.vars=c("IID", "deathDate"), measure.vars=patterns("specimenDate_*", "SARS_CoV2_"),
              variable.name="testNumber", value.name=c("testDate", "testResult"))
covid <- covid[!is.na(IID) & !is.na(testDate) & testDate != ""]
covid[, testDate := as.IDate(testDate, format="%d%b%Y")]
covid[, testNumber := as.integer(as.vector(testNumber))]
covid[, deathDate := as.IDate(deathDate, format="%d%b%Y")]

# Collate information about the sometimes multiple tests per person
collated_covid <- covid[, .(n_covid_tests = .N,
  first_covid_test = min(testDate),
  last_covid_test = max(testDate),
  any_positive = any(testResult == 1)
), by = .(IID)]

# Get date of first positive covid test
positive <- covid[testResult == 1]
first_positive <- positive[,.SD[which.min(testDate)], by=IID]
collated_covid[first_positive, on = .(IID), first_positive := i.testDate]

# Get date of laste positive covid test
last_positive <- positive[,.SD[which.max(testDate)], by=IID]
collated_covid[last_positive, on = .(IID), last_positive := i.testDate]

# Get date of first negative test after all positive tests
negative_after_positive <- covid[last_positive, on = .(IID, testDate > testDate), 
  .(IID, last_positive=i.testDate, testDate=x.testDate, testNumber), nomatch=0]
first_negative <- negative_after_positive[,.SD[which.min(testDate)], by=IID]
collated_covid[first_negative, on = .(IID), first_negative := i.testDate]

# Get time between first positive test and first negative test post covid.
# Note most of these are people who tested negative several months after the
# first peak - not negative tests following first peak covid infection
collated_covid[, covid_duration := first_negative - first_positive]

# Get deaths following positive covid results
covid_deaths <- covid[IID %in% positive$IID & !is.na(deathDate)]
collated_covid[covid_deaths, on = .(IID), death_after_covid := deathDate]

# Deaths within 28 days of last positive test result
collated_covid[(death_after_covid - last_positive) <= 28, death_within_28d := death_after_covid]

# record in HES whether event is pre or post covid
pandemic_start <- collated_covid[,min(first_covid_test)]
hes[, pandemic_before_event := ifelse(eventDate >= pandemic_start, 1L, 0L)]

# Record whether event is pre or post covid infection
hes[, covid_before_event := 0L]
hes[collated_covid[(any_positive)], on = .(IID), covid_before_event := ifelse(eventDate >= first_positive, 1L, 0L)]

# Extract loaded results out of local script environment
}) # end environment execution

hes <- hes_script_env$hes
covid <- hes_script_env$collated_covid
rm(hes_script_env)

