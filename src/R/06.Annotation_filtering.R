######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: Annotation cleaning and export data.
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################
library(tidymass)
library(tidyverse)

message(crayon::green("load annotated objects. please wait..."))
load("object_neg_anno.rds")
load("object_pos_anno.rds")
load("~/.HPC_tidymass/MS_db/labID2INCHIKEY.rda")
message(crayon::green("Running original filtering!"))
# original annotation output ----------------------------------------------
dir.create("Original_annotation",showWarnings = F,recursive = T)
object_merge_original <- 
  merge_mass_dataset(
    x = object_neg_anno,
    y = object_pos_anno,
    sample_direction = "inner",
    variable_direction = 'full',
    sample_by = 'sample_id',
    variable_by = c("variable_id","mz","rt")
  ) %>% 
  activate_mass_dataset("sample_info") %>% 
  filter(class != "QC")

ori_vari_info <-
  object_merge_original %>% 
  extract_variable_info() %>% 
  select(variable_id,mz,rt)

anno_original_tbl <-
  object_merge_original %>% 
  extract_annotation_table() %>% 
  group_by(variable_id,Compound.name) %>% 
  slice_max(order_by = Total.score,n =1 ) %>% 
  left_join(.,lab2inchi,by = "Lab.ID") %>% 
  left_join(.,ori_vari_info,by = "variable_id")


ori_exp_mat <- 
  object_merge_original %>% 
  extract_expression_data() %>% 
  rownames_to_column("Compound_ID")
ori_sample_info <-
  object_merge_original %>% 
  extract_sample_info() 
ori.out = list(
  Full = anno_original_tbl,
  Exp_mat = ori_exp_mat,
  sample_info = ori_sample_info
)
writexl::write_xlsx(x = ori.out,path = "Original_annotation/Original_metadata.xlsx")
save(object_merge_original,file = "Original_annotation/object_merge_original.rds")


# filter out --------------------------------------------------------------
message(crayon::green("Running rough filtering!"))
dir.create("Clean_annotation",showWarnings = F,recursive = T)
dir.create("Clean_annotation/With_MS1")

object_filter_rough <-
  object_merge_original %>% 
  activate_mass_dataset("annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  group_by(Compound.name) %>% 
  filter(Level == min(Level)) %>% 
  filter(Total.score == max(Total.score)) %>% 
  group_by(variable_id) %>% 
  filter(Total.score == max(Total.score)) %>% 
  slice_head(n = 1) %>% 
  filter(Adduct == "(M-H)-" |
           Adduct == "(M-H2O-H)-"|
           Adduct == "(M+NH4-2H)-"|
           Adduct == "(2M-H)-"|
           Adduct == "(M+F)-" |
           Adduct == "(M+H)+" |
           Adduct == "(2M+H)+"|
           Adduct == "(M+H-H2O)+"|
           Adduct == "(M+NH4)+"|
           Adduct == "(M+H-2H2O)+" |
           Adduct == "(M+H-2H2O)+"|
           Adduct == "(2M+NH4)+"
           )

rough_vari_info <-
  object_filter_rough %>% 
  extract_variable_info() %>% 
  select(variable_id,mz,rt)

rough_exp_mat <- 
  object_filter_rough %>% 
  extract_expression_data() %>% 
  rownames_to_column("Compound_id")

rough_sample_info <-
  object_filter_rough %>% 
  extract_sample_info()

anno_rough_filter = 
  object_filter_rough %>% 
  extract_annotation_table() %>% 
  unique() %>% 
  left_join(.,lab2inchi,by = "Lab.ID") %>% 
  left_join(.,rough_vari_info,by = "variable_id")

out.rough <- list(
  exp_mat = rough_exp_mat,
  anno = anno_rough_filter,
  sample_info = rough_sample_info
)

writexl::write_xlsx(x = ori.out,path = "Clean_annotation/With_MS1/Original_metadata.xlsx")
save(object_filter_rough,file = "Clean_annotation/With_MS1/object_filter_rough.rds")

##> remove redundant features.
message(crayon::green("Remove redundant filtering!"))
dir.create("Clean_annotation/rm_redundant",showWarnings = F,recursive = T)

object_filter_rm_redundant <-
  object_filter_rough %>% 
  activate_mass_dataset("annotation_table") %>% 
  group_by(Compound.name_case = stringi::stri_trans_totitle(Compound.name)) %>% 
  slice_head(n = 1)

rmd_vari_info <-
  object_filter_rm_redundant %>% 
  extract_variable_info() %>% 
  select(variable_id,mz,rt)

rmd_exp_mat <- 
  object_filter_rm_redundant %>% 
  extract_expression_data() %>% 
  rownames_to_column("Compound_id")

rmd_sample_info <-
  object_filter_rm_redundant %>% 
  extract_sample_info()

anno_rmd_filter = 
  object_filter_rm_redundant %>% 
  extract_annotation_table() %>% 
  unique() %>% 
  left_join(.,lab2inchi,by = "Lab.ID") %>% 
  left_join(.,rmd_vari_info,by = "variable_id")

out.rmd <- list(
  exp_mat = rmd_exp_mat,
  anno = anno_rmd_filter,
  sample_info = rmd_sample_info
)

writexl::write_xlsx(x = out.rmd,path = "Clean_annotation/rm_redundant/rm_redundant_metadata.xlsx")
save(object_filter_rm_redundant,file = "Clean_annotation/rm_redundant/object_filter_rmd.rds")

##> Keep MS2 only

message(crayon::green("Select hight confidance compounds!"))
dir.create("Clean_annotation/Only_MS2",showWarnings = F,recursive = T)

object_filter_high_confidence <-
  object_filter_rm_redundant %>% 
  activate_mass_dataset("annotation_table") %>% 
  filter(Level != 3)

high_vari_info <-
  object_filter_high_confidence %>% 
  extract_variable_info() %>% 
  select(variable_id,mz,rt)

high_exp_mat <- 
  object_filter_high_confidence %>% 
  extract_expression_data() %>% 
  rownames_to_column("Compound_id")

high_sample_info <-
  object_filter_high_confidence %>% 
  extract_sample_info()

anno_high_filter = 
  object_filter_high_confidence %>% 
  extract_annotation_table() %>% 
  unique() %>% 
  left_join(.,lab2inchi,by = "Lab.ID") %>% 
  left_join(.,high_vari_info,by = "variable_id")

out.high <- list(
  exp_mat = high_exp_mat,
  anno = anno_high_filter,
  sample_info = high_sample_info
)

writexl::write_xlsx(x = out.high,path = "Clean_annotation/Only_MS2/high_confidence_metadata.xlsx")
save(object_filter_high_confidence,file = "Clean_annotation/Only_MS2/object_filter_high.rds")


