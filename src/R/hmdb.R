#################################
#       Prj: Tidymass
#       Assignment: Classyfire
#       Author: shawn
#       Date: Jun 1, 2022
#       Location: HENU
##################################

library(tidyverse)
library(MDAtoolkits)
library(webchem)
library(classyfireR)

load("/home/data/public/01.database/01.Database/02.meta_db/hmdb_database0.0.3.rda")



hmdb_compound_db <- 
  hmdb_database0.0.3@spectra.info %>% 
  select(Lab.ID,Compound.name,SMILES,mz,Formula,RT) 

Ident_mat = hmdb_compound_db %>% 
  select(Lab.ID,SMILES,Formula,mz,RT) %>% 
  mutate(Judge = T) %>% 
  setNames(c("Compound_ID","name","mf","mw","RT","Judge"))
step1_result= mda_data_org(compound_info = Ident_mat,source = "other")
head(step1_result)
step2_result = mda_get_cid(
  data_info = step1_result,
  type = 'multiple',match = 'all',core_num = 50,from = 'smile'
)

step3_result = mda_pubchem_crawler(
  cid_info = step2_result,type = 'multiple',core_num = 50
)

step4_result = mda_classfire(
  query = step3_result,type = "multiple"
)

step5_result = mda_CTS_kegg(
  query = step3_result,type = "multiple"
)
colnames(step1_result)[1] = "compound_id"
step1_result %>% 
  mutate(mf = gsub(" ","",mf)) -> step1_result
step6_result = mda_merge_info(
  orz_data = step1_result,name2cid = step2_result,pubchem_detail = step3_result,pubchem_class = step4_result,CTS_kegg = step5_result
)


step7_result = mat_all_final_merge %>% 
  select(
    Compound_ID,Name,val_mean,level,level2,ions,Checked,Formula,`Molecular Weight`,`RT [min]`
  ) %>% 
  rename("compound_id" = "Compound_ID") %>% 
  inner_join(step6_result,by = "compound_id") %>% 
  rename("Compound_name" = "query") %>% 
  relocate(Compound_name,.after = compound_id) %>% 
  rename("CD_Check" = "Checked") 
colnames(step7_result)
step8_result = step7_result %>% 
  filter(High_identical == "TRUE")
save.image(file = "img.Rdata")
write.csv(step7_result,"01.Metbo/Table1_Compound_Ident_Class.csv",row.names = F)