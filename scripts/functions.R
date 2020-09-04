# function to map WormJam InChI key to external DB id ==========================
wormJam_mapper <- function(inchikey, mapper) {
  
  # get all codes required in WormJam
  wormjam_codes <- .get_wormjam_codes()
  
  # generate empty list for ids
  ids <- list()
  
  for(wormjam_code in names(wormjam_codes)) {
    
    wormjam_mapping <- BridgeDbR::map(mapper,
                                      inchikey,
                                      source = getSystemCode("InChIKey"),
                                      target = wormjam_codes[[wormjam_code]])
    
    if(length(wormjam_mapping) > 0) {
      
      # isolate based on regex
      if(wormjam_code == "ChEBI") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^CHEBI:\\d+$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "KEGG") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^C\\d+$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "MetaCyc") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^CPD-\\d{5}$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "HMDB") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^HMDB\\d+$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "LipidMaps") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^LM(FA|GL|GP|SP|ST|PR|SL|PK)[0-9]{4}([0-9a-zA-Z]{4,6})?$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "SwissLipids") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^SLM:\\d+$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "Wikidata") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^Q\\d+$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "PubChem") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^\\d+$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "Metabolights") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^MTBLC\\d+$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } else if(wormjam_code == "Chemspider") {
        
        ids_list <- unlist(str_extract_all(wormjam_mapping$mapping,
                                           "^\\d+$"))
        
        ids[[wormjam_code]] <- paste0(ids_list, collapse = ";")
        
      } 
      
    } else {
      
      ids[[wormjam_code]] <- NA
      
    }
  }
  
  return(ids)
  
}

# helper function for BridgDbR codes in WormJam ================================
.get_wormjam_codes <- function() {
  wormjam_codes <- list(
    "ChEBI" = BridgeDbR::getSystemCode("ChEBI"),
    "KEGG" = BridgeDbR::getSystemCode("KEGG Compound"),
    "MetaCyc" = BridgeDbR::getSystemCode("MetaCyc"),
    "HMDB" = BridgeDbR::getSystemCode("HMDB"),
    "LipidMaps" = BridgeDbR::getSystemCode("LIPID MAPS"),
    "SwissLipids" = BridgeDbR::getSystemCode("SwissLipids"),
    "Wikidata" = BridgeDbR::getSystemCode("Wikidata"),
    "PubChem" = BridgeDbR::getSystemCode("PubChem-compound"),
    "Metabolights" = BridgeDbR::getSystemCode("MetaboLights Compounds"),
    "Chemspider" = BridgeDbR::getSystemCode("Chemspider")
  )
  
  return(wormjam_codes)
}  
