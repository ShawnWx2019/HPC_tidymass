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
message(crayon::green("Running original filtering!"))

# annotation filtering ----------------------------------------------------
##> variable data
##> 
object_neg_anno.respect = object_neg_list[[1]]
object_pos_anno.respect = object_pos_list[[1]]
vari_info_neg <-
  object_neg_anno.respect %>% 
  extract_variable_info() %>% 
  select(variable_id,mz,rt,contains("na_freq"),norm.rsd)
vari_info_pos <-
  object_pos_anno.respect %>% 
  extract_variable_info() %>% 
  select(variable_id,mz,rt,contains("na_freq"),norm.rsd)
vari_info <- rbind(vari_info_pos,vari_info_neg)
##> data filtering
feature_annotation_all <- 
  purrr::map_dfr(.x = 1:length(object_pos_list),.f = function(.x){
    anno_tbl.pos = object_pos_list[[.x]] %>% extract_annotation_table() %>% dplyr::select(-ms2_files_id) %>% 
      filter(Adduct == "(M+H)+")
    anno_tbl.neg = object_neg_list[[.x]] %>% extract_annotation_table() %>% dplyr::select(-ms2_files_id) %>% 
      filter(Adduct == "(M-H)-")
    out = rbind(anno_tbl.pos,anno_tbl.neg)
    return(out)
  }) %>% distinct() %>% arrange(variable_id)
feature_annotation_all <-
  left_join(vari_info,feature_annotation_all,by = "variable_id")

dir.create("Original_annotation",showWarnings = F,recursive = T)

writexl::write_xlsx(x = feature_annotation_all,path = "Original_annotation/01.Original_FeatureAnnotation.xlsx")

expmat_all_clean_feature.pos <- 
  object_pos_anno.respect %>% 
    extract_expression_data() %>% 
    rownames_to_column("variable_id")

expmat_all_clean_feature.neg <- 
  object_neg_anno.respect %>% 
  extract_expression_data() %>% 
  rownames_to_column("variable_id")
if(ncol(expmat_all_clean_feature.pos) > ncol(expmat_all_clean_feature.neg)) {
  expmat_all_clean_feature.pos = 
    expmat_all_clean_feature.pos %>% select(colnames(expmat_all_clean_feature.neg))
} else if(
  ncol(expmat_all_clean_feature.pos) < ncol(expmat_all_clean_feature.neg)
) {
  expmat_all_clean_feature.neg = 
    expmat_all_clean_feature.neg %>% select(colnames(expmat_all_clean_feature.pos))
}
expmat_all <- rbind(expmat_all_clean_feature.pos,expmat_all_clean_feature.neg)

writexl::write_xlsx(x = expmat_all,path = "Original_annotation/02.Original_FeatureAccumulationMatrix.xlsx")
                    
Sample_info_all_clean_feature <- 
  object_neg_anno.respect %>% 
  extract_sample_info() 

writexl::write_xlsx(x = Sample_info_all_clean_feature,path = "Original_annotation/03.Original_Sample_info.xlsx")


# Represent feature filtring ----------------------------------------------

Represent_feature_anno <- 
  feature_annotation_all %>% 
  group_by(variable_id) %>% 
  filter(Level == min(Level)) %>% 
  slice_max(order_by = Total.score,n = 1) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  group_by(Compound.name) %>% 
  slice_max(order_by = Total.score,n = 1) %>% 
  slice_head(n = 1) %>% 
  arrange(variable_id)

dir.create("Feature_filtering",showWarnings = F,recursive = T)

writexl::write_xlsx(x = Represent_feature_anno,path = "Feature_filtering/01.Filtered_FeatureAnnotation.xlsx")

expmat_filter <- 
  inner_join(
   Represent_feature_anno %>% select(variable_id),expmat_all
  )

writexl::write_xlsx(x = expmat_filter,path = "Feature_filtering/02.Filtered_FeatureAccumulationMatrix.xlsx")
writexl::write_xlsx(x = Sample_info_all_clean_feature,path = "Feature_filtering/03.Filtered_Sample_info.xlsx")

