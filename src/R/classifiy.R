#################################
#       Prj: Tidymass
#       Assignment: Classyfire
#       Author: shawn
#       Date: Jun 1, 2022
#       Location: HENU
##################################

load("~/SynologyDrive/database/02.MS/MSdb/labID2INCHIKEY.rda")
library(tidyverse)
library(MDAtoolkits)
lab2inchi_ok <-
  lab2inchi %>% 
  filter(!is.na(INCHIKEY))
inchikey_uniq <- 
  lab2inchi_ok %>% 
  select(INCHIKEY) %>% unique() %>% 
  setNames("InChIKey")

write.table(x = inchikey_uniq,file = "demo/inchikey_alldb.xls",row.names = F,sep = "\t")




load("~/SynologyDrive/database/02.MS/MSdb/fiehn_hilic_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/hmdb_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/kegg_ms1_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/massbank_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/mona_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/snyder_database_hilic0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/KNApSAcK_ms1_database.rda")
load("~/SynologyDrive/database/02.MS/MSdb/plantcyc_ms1_database0.0.2.rda")
load("~/SynologyDrive/database/02.MS/MSdb/RIKEN_PlaSMA_database0.0.1.rda")
load("~/SynologyDrive/database/02.MS/MSdb/Natural_products_database_v1.rds")


libID2INCHIKEY = function(db){
  tryCatch({
    out =
      db@spectra.info %>%
      select(Lab.ID,
             matches('compound.name',ignore.case = T),
             matches("pubchem.id",ignore.case = T),
             matches("inchikey",ignore.case = T),
             matches('formula',ignore.case = T),
             mz,
             matches('cas.id',ignore.case = T)
             ) %>%
      distinct()
    if(ncol(out) == 5) {
      final = out %>% set_names("Lab.ID","Compound.name","PubChem.ID","INCHIKEY","Formula","mz","CAS.ID")
    } else if (
      isTRUE(
        str_detect(
          string = colnames(out)[3],
          pattern = regex(pattern = "pubchem",ignore_case = T)
        )
      )) {
      final = data.frame(
        Lab.ID = out[1],
        Compound.name = out[2],
        PubChem.ID = out[3],
        INCHIKEY = NA,
        Formula = out[4],
        mz = out[5],
        CAS.ID = out[6]
      )
    } else if (
      isTRUE(
        str_detect(
          string = colnames(out)[3],
          pattern = regex(pattern = "inchikey",ignore_case = T)
        )
      )) {
      final = data.frame(
        Lab.ID = out[1],
        Compound.name = out[2],
        PubChem.ID = NA,
        INCHIKEY = out[3],
        Formula = out[4],
        mz = out[5],
        CAS.ID = out[6]
      )
    }
  })
  colnames(final) = c("Lab.ID","Compound.name","PubChem.ID","INCHIKEY","Formula","mz","CAS.ID")
  final$INCHIKEY = gsub(pattern = "InChIKey=",replacement = "", x = final$INCHIKEY,ignore.case = T)
  final$Formula = gsub(pattern = " ","",final$Formula)
  return(final)
}

tmp1 = libID2INCHIKEY(db = KNApSAcK_ms1_database0.0.1)
tmp2 = libID2INCHIKEY(db = fiehn_hilic_database0.0.3)
tmp3 = libID2INCHIKEY(db = mona_database0.0.3)
tmp4 = libID2INCHIKEY(db = hmdb_database0.0.3)
tmp5 = libID2INCHIKEY(db = kegg_ms1_database0.0.3)
tmp6 = libID2INCHIKEY(db = massbank_database0.0.3)
tmp7 = libID2INCHIKEY(db = RIKEN_PlaSMA_database0.0.1)
tmp8 = libID2INCHIKEY(db = plantcyc_ms1_database0.0.2)
tmp9 = libID2INCHIKEY(db = Natural_products_database_v1)

lab2inchi =
  rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9) %>% as_tibble()
