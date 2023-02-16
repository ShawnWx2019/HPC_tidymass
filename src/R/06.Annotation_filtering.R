######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: Annotation cleaning and export data.
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################
library(tidymass)
library(tidyverse)

setwd("demo/Annotation/")
load("object_neg_anno.rds")
load("object_pos_anno.rds")


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

anno_original_tbl <-
  object_merge_original %>% 
  extract_annotation_table() %>% 
  group_by(variable_id,Compound.name) %>% 
  slice_max(order_by = Total.score,n =1 )

anno_original_tbl_lv2 <-
  anno_original_tbl %>% 
  filter(Level != 3)
ori_exp_mat <- 
  object_merge_original %>% 
  extract_expression_data() %>% 
  rownames_to_column("Compound_ID")
object_pos_anno %>% extract_sample_info()
ori_sample_info <-
  object_merge_original %>% 
  extract_sample_info() 
ori.out = list(
  Full = anno_original_tbl,
  level2 = anno_original_tbl_lv2,
  Exp_mat = ori_exp_mat,
  sample_info = ori_sample_info
)
writexl::write_xlsx(x = ori.out,path = "Original_annotation/Original_metadata.xlsx")
save(object_merge_original,file = "Original_annotation/object_merge_original.rds")


# filter out --------------------------------------------------------------

dir.create("Clean_annotation",showWarnings = F,recursive = T)

