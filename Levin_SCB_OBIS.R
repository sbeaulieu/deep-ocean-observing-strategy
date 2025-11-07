#Levin_SCB_OBIS
#Stace Beaulieu
#2025-11-06

# R script to standardize Lisa Levin's SCB data to Darwin Core (DwC)
# and output tables for OBIS
#   event Core table
#   occurrence extension table
#   (and consider emof table)

#Setting wd

setwd("/Users/sbeaulieu/Downloads")

#Load required packages

library(readxl)
library(dplyr)
library(tidyr)
library(worrms) # use worrms for QC and to populate taxonRank
library(stringr) # when using worrms for QC

#Load datasheets

counts_input <- readxl::read_xlsx("levin_obis_data_SCB Hardground macrofauna.xlsx", sheet = "Occurrence Info")
#counts_input <- readxl::read_xlsx("Levin_Guraieb SCB Hard Substrates_ SCB.xlsx", sheet = "Macrofauna Counts")
# 2 more than WoRMS table bc manually removed

events_input <- readxl::read_xlsx("levin_obis_data_SCB Hardground macrofauna.xlsx", sheet = "eventTable Info")
#events_input <- readxl::read_xlsx("Levin_Guraieb SCB Hard Substrates_ SCB.xlsx", sheet = "Depth LL  T O2 Substrate etc.")

# taxa sheet edited from WoRMS Taxon Match tool output
taxa_input <- readxl::read_xlsx("levin_obis_data_SCB Hardground macrofauna.xlsx", sheet = "WoRMS match")

# some edits are needed for the taxa_input
# I want to use the column LSID for scientificNameID
# note some LSID entries need edits
# need to strip LSID prefix for wm_id2name
stripped <- sub(".*:(\\d+)$", "\\1", taxa_input$LSID)
stripped_int <- as.integer(stripped[grepl("^[0-9]+$", stripped)])
check_names <- wm_id2name_(id = stripped_int)
#need to edit the LSIDs input first so list same length
taxa_input_filtered <- taxa_input %>%
  filter(str_starts(LSID, "urn"))
taxa_input_filtered$check_names <- unlist(check_names) 
taxa_input_filtered$equal <- taxa_input_filtered$ScientificName == taxa_input_filtered$check_names
taxa_input_filtered[!taxa_input_filtered$equal, ]


#Initiate DwC events table
# rename columns to DwC terms
event_dwc <- events_input %>%
  rename(verbatimLabel = "Sample number") %>%
  rename(decimalLatitude = "Latitude") %>%
  rename(decimalLongitude = "Longitude") %>%
  rename(habitat = "Habitat") %>%
  rename(verbatimDepth = "Water Depth (M)") %>%
  rename(sampleSizeValue = "Rock Surface Area (cm2)") # note one value is character E

# add DwC sampleSizeUnit
event_dwc$sampleSizeUnit <- "Rock Surface Area (cm2)"

# create eventID by concatenating cruise, dive, rock, number, and "Sample number"
# "Sample number" will be redundant for some but not all cruises
event_dwc <- event_dwc %>%
  unite(col='eventID', c('Cruise', 'Dive', 'Rock number', 'verbatimLabel'), sep = "_", remove = FALSE, na.rm = FALSE)

# create DwC eventRemarks from multiple columns
# include temp and oxygen but ultimately move to emof
# note this could be improved by adding descriptive strings for each column's value
event_dwc <- event_dwc %>%
  unite(col='eventRemarks', c('Station', 'Temperature', 'Oxygen (umol/L)'), sep = "|", remove = TRUE, na.rm = FALSE)

# create min and max Depth equiv to verbatimDepth
event_dwc$minimumDepthInMeters <- event_dwc$verbatimDepth
event_dwc$maximumDepthInMeters <- event_dwc$verbatimDepth

# consider adding countryCode needed for GBIF
event_dwc$countryCode <- "US"

#Finish DwC events table
# remove extra columns
# DwC basisOfRecord will go into occurrence table
# consider publishing Temperature and Oxygen in emof table
event <- select(event_dwc, -basisOfRecord, -Cruise, -Dive, -'Rock number', - 'Proximity to shore',-starts_with(".."))
#write.csv(event, 'event.csv', row.names=FALSE)
               
# consider adding geodeticDatum

#Initiate DwC occurrence extension table
# first add column with row counter for "verbatimID" (was "Species morphotype")
# to be used as suffix for occurrenceID
counts <- counts_input %>% mutate(vIrow = row_number())
counts$vIrow <- sprintf("%03d", counts$vIrow) # pad with leading zero to 3 characters

# QC that the counts column headers
# are equivalent to the verbatimLabel 
colnames_counts <- colnames(counts)
IDs_event <- event$verbatimLabel
matches <- colnames_counts %in% IDs_event
matches
# the first and last are FALSE bc verbatimID and vIrow

# transform wide to long
counts_long <- counts %>%
  pivot_longer(
    cols = 2:83, # new column vIrow was added to far right
    names_to = "verbatimLabel",
    values_to = "individualCount"
  )

# consider removing zeroes at this point in script
# retain occurrenceStatus present
counts_long_present <- filter(counts_long, individualCount > 0)
# QC with sum
sum(counts_long_present$individualCount)
sum(counts_long$individualCount)

# add eventID into occurrence table
# but remove extra event columns
counts_long_present_eventID <- full_join(counts_long_present, event, by = 'verbatimLabel')
counts_long_present_eventID <- counts_long_present_eventID[, 1:5]

# rename column and add DwC term
occurrence_dwc <- counts_long_present_eventID %>%
  rename(verbatimIdentification = "verbatimID") %>%
  unite(occurrenceID, c("eventID", "vIrow"), sep = "_", remove = FALSE)

# join taxonomic standardization

taxa <- taxa_input %>%
  rename(verbatimIdentification = "verbatimID")
occurrence_dwc <- full_join(occurrence_dwc, taxa, by = 'verbatimIdentification')

# add remaining required DwC terms for OBIS
occurrence_dwc$basisOfRecord <- "PreservedSpecimen"

# already removed zero individualCount
occurrence_dwc$occurrenceStatus <- "present" # if remove zero individualCount

# rename to DwC and exclude non-DwC columns
occurrence_dwc <- occurrence_dwc %>%
  rename(scientificName = "ScientificName") %>%
  rename(scientificNameID = "LSID") %>%
  rename(kingdom = "Kingdom")

col_order <- c("occurrenceID", "verbatimIdentification", "individualCount", "scientificName", "scientificNameID",
               "kingdom", "occurrenceStatus",
               "basisOfRecord", "eventID")
occurrence <- occurrence_dwc[, col_order]
#write.csv(occurrence, 'occurrence.csv', row.names=FALSE)

# QC matching bw event and occurrence
occ_event_match <- full_join(occurrence, event, by = 'eventID')