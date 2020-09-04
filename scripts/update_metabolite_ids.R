# load required libraries ------------------------------------------------------
library(BridgeDbR)
library(tidyverse)

# load functions required for this script --------------------------------------
source("scripts/functions.R")

# general setup required for script --------------------------------------------
# load the current metabolomics BridgeDB
mapper <- BridgeDbR::loadDatabase("E:/Databases/bridgedb/2020-08-09/metabolites_20200809.bridge")

# list all files in literature curation ----------------------------------------
literature_files <- list.files("literature",
                               pattern = ".txt$",
                               full.names = TRUE)

literature_metabolites <- tibble()

# iterate over literature files
for(literature_file in literature_files) {
  
  metabolites <- read_tsv(literature_file)
  
  # iterate through all compounds in table
  for(i in 1:nrow(metabolites)) {
    
    # get inchikey for mapping to IDs
    inchikey <- metabolites$InChIKey[i]
    
    if(!is.na(inchikey)) {
      
      # perform mapping
      ids <- wormJam_mapper(inchikey, mapper)
      
      print(ids)
      
      # add ids to compound table
      metabolites$ChEBI[i] <- ids[["ChEBI"]]
      metabolites$KEGG[i] <- ids[["KEGG"]]
      metabolites$BioCyc[i] <- ids[["MetaCyc"]]
      metabolites$HMDB[i] <- ids[["HMDB"]]
      metabolites$LipidMaps[i] <- ids[["LipidMaps"]]
      metabolites$SwissLipids[i] <- ids[["SwissLipids"]]
      metabolites$Wikidata[i] <- ids[["Wikidata"]]
      metabolites$PubChem[i] <- ids[["PubChem"]]
      metabolites$Metabolights[i] <- ids[["Metabolights"]]
      metabolites$ChemSpider[i] <- ids[["Chemspider"]]
    }
  }
  
  # sort out only unique rows and add to literature list
  metabolites <- distinct(metabolites)
  literature_metabolites <- bind_rows(literature_metabolites, metabolites)
  
  # write result file
  write_tsv(metabolites, literature_file)
  
}

# filter out only unique values for literature metabolites
literature_metabolites_unique <- distinct(literature_metabolites)
