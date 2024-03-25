####################################################################################
#           Prj: HPC-tidymass
#           Assignment: mv redun between pos and neg
#           Author: Shawn Wang
#           Date: Wed Jan 03, 2024
####################################################################################
suppressMessages(library(tidymass))
suppressMessages(library(tidyverse))
suppressMessages(library(furrr))
load("object_neg.ms2.rds")
load("object_pos.ms2.rds")
suppressMessages(library(crayon))
options(warn = -1)
msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic
# compare_features_ms2 = function(object.pos = object_pos.ms2,
#          object.neg = object_neg.ms2,
#          feature_pos = "M209T132_POS",
#          feature_neg = "M207T131_NEG",
#          n_peak=20) {
#   ms2_pos_data <- object_pos.ms2 %>% extract_ms2_data() %>% setNames('ms2_data')
#   ms2_spec.pos = ms2_pos_data$ms2_data@ms2_spectra
#   index.pos = data.frame(
#     variable = ms2_pos_data$ms2_data@variable_id
#   ) %>% 
#     mutate(order = 1:nrow(.))
#   ms2_neg_data <- object_neg.ms2 %>% extract_ms2_data() %>% setNames('ms2_data')
#   ms2_spec.neg = ms2_neg_data$ms2_data@ms2_spectra
#   index.pos = data.frame(
#     variable = ms2_pos_data$ms2_data@variable_id
#   ) %>% 
#     mutate(order = 1:nrow(.))
#   index.neg = data.frame(
#     variable = ms2_neg_data$ms2_data@variable_id
#   ) %>% 
#     mutate(order = 1:nrow(.))
#   check_index = index.pos %>% 
#     filter(variable == feature_pos) %>% pull(order) %>% as.numeric()
#   check_index2 = index.neg %>% 
#     filter(variable == feature_neg) %>% pull(order) %>% as.numeric()
  
#   check1 = ms2_spec.pos[[check_index]] %>%
#     as.data.frame %>% mutate(intensity = intensity/max(intensity)*100,
#                              mz = mz - 1) %>% 
#     slice_max (order_by = intensity,n = n_peak)
#   check2 = ms2_spec.neg[[check_index2]] %>%
#     as.data.frame %>% mutate(intensity = intensity/max(intensity)*100,
#                              mz = mz + 1) %>% 
#     slice_max (order_by = intensity,n = n_peak)
#   check1_df <- 
#     check1 %>% 
#     mutate(
#       variable_id = feature_pos,
#       x = mz,
#       xend = mz,
#       y = 0,
#       yend = intensity,
#       tag = case_when(
#         intensity > 10 ~ mz,
#         TRUE ~ NA
#       )
#     )
#   check2_df <- 
#     check2 %>% 
#     mutate(
#       variable_id = feature_neg ,
#       x = mz,
#       xend = mz,
#       y = 0,
#       yend = -intensity,
#       tag = case_when(
#         intensity > 10 ~ mz,
#         TRUE ~ NA
#       )
#     )
  
#   check_merge = rbind(check1_df,check2_df)
  
#   p = ggplot(check_merge,
#              mapping = aes(x = x,xend = xend,y = y,yend = yend,color = variable_id)) + 
#     geom_segment(size = 1)+
#     geom_text(mapping = aes(x = x,y = yend,label = tag),color = 'black')+
#     theme_bw() + 
#     scale_y_continuous(breaks = c(-100, -50, 0, 50, 100),labels = c('100',"50","0","50","100")) +
#     xlab("fixed mz")+
#     ylab('intensity') +
#     theme(
#       axis.ticks = element_line(size = 1),
#       panel.border = element_rect(linewidth = 1),
#       axis.text = element_text(size = 12)
#     )
#   return(p)
# }
##> 
pos_final_rest <- readxl::read_xlsx("RemoveRedundancy/POS_final_rest.xlsx")
neg_final_rest <- readxl::read_xlsx("RemoveRedundancy/NEG_final_rest.xlsx")

future::plan('multicore',workers = 16)
message(msg_run('group features my mz_tol = 1 da'))
similar_posvsneg <- 
  furrr::future_map_dfr(.x = 1:nrow(pos_final_rest),.f = function(.x){
    mz_pos = pos_final_rest[.x,4] %>% as.numeric()
    mz_pos_fix = mz_pos - 1
    rt_pos = pos_final_rest[.x,5]%>% as.numeric()
    mz_range = c(mz_pos_fix - 1,mz_pos_fix + 1)
    
    tmp_neg_final_rest <-
      neg_final_rest %>% 
      dplyr::filter((mz.ms2 > mz_range[1]-1) & (mz.ms2 < mz_range[2]-1))
    
    if(nrow(tmp_neg_final_rest) > 0) {
      out =  map_dfr(1:nrow(tmp_neg_final_rest),.f = function(.y) {
        data.frame(
          Pos.ID = pos_final_rest[.x,1] %>% as.character(),
          Neg.ID = tmp_neg_final_rest[.y,1] %>% as.character(),
          mz_pos = mz_pos,
          mz_neg = tmp_neg_final_rest[.y,4] %>% as.numeric(),
          mz_diff = abs(mz_pos_fix - tmp_neg_final_rest[.y,4] %>% as.numeric() + 1) * (10^6) / mz_pos_fix,
          rt_pos = rt_pos,
          rt_neg = tmp_neg_final_rest[.y,5] %>% as.numeric(),
          rt_diff = abs(rt_pos - tmp_neg_final_rest[.y,5] %>% as.numeric())
        )
      })
    } else {
      out = data.frame(
        Pos.ID = NA,
        Neg.ID = NA,
        mz_pos = NA,
        mz_neg = NA,
        mz_diff = NA,
        rt_pos = NA,
        rt_neg = NA,
        rt_diff = NA
      )
    }
  },.progress = T) %>% 
    drop_na()



# compare_features_ms2(feature_pos = "M113T63_2_POS",object.pos = "M111T63_NEG",n_peak = 20)
##> MS2 没有可比性 比较没有意义
message(msg_run('Select the pick with the highst intensity.'))
redundant_group <-
  similar_posvsneg %>% 
  tibble() %>% 
    filter(mz_diff < 100000 & rt_diff < 3) %>% 
    group_by(round(mz_pos,0),round(mz_neg,0)) %>% 
    left_join(pos_final_rest %>% select(variable_id,mean) %>% setNames(c('Pos.ID','Pos.mean'))) %>% 
    left_join(neg_final_rest %>% select(variable_id,mean) %>% setNames(c('Neg.ID','Neg.mean'))) %>% 
    filter((Pos.mean + Neg.mean) == max(Pos.mean + Neg.mean)) %>% 
    mutate(
      losser = case_when(
        Pos.mean > Neg.mean ~ Neg.ID,
        TRUE ~ Pos.ID
      )
    ) %>% 
    ungroup() %>% 
    relocate(losser,.after = Neg.ID)

neg_final <-
neg_final_rest %>% 
anti_join(
  data.frame(
    variable_id = redundant_group %>% filter(str_detect(losser,'NEG')) %>% pull(losser)
  )
)

pos_final <- 
pos_final_rest %>% 
anti_join(
  data.frame(
    variable_id = redundant_group %>% filter(str_detect(losser,'POS')) %>% pull(losser)
  )
)

feature_final <- 
  rbind(pos_final,neg_final)
Sys.sleep(3)
message(msg_run('export non-redundant features'))
writexl::write_xlsx(
  list(
    negative = neg_final,
    positive = pos_final,
    all = feature_final
  ),"RemoveRedundancy/Non_redundant_feature.xlsx"
)

## export annotation    =====================================

## load annotation project =====================================
Sys.sleep(3)

message(msg_run('export metadata'))
load("object_neg_anno.rds")
load("object_pos_anno.rds")

pos_variables <- pos_final
neg_variables <- neg_final

## export expression matrix   =====================================
object_clean.pos <- 
    object_pos.ms2 %>% 
    activate_mass_dataset('variable_info') %>% 
    left_join(pos_variables %>% select(-mz,-rt),by = 'variable_id') %>% 
    filter(tag == 'marked')

object_clean.neg <- 
    object_neg.ms2 %>% 
    activate_mass_dataset('variable_info') %>% 
    left_join(neg_variables %>% select(-mz,-rt),by = 'variable_id') %>% 
    filter(tag == 'marked')

object_clean <- merge_mass_dataset(
    x = object_clean.pos,y = object_clean.neg,sample_direction = 'inner', variable_direction = 'full',variable_by = c(
        'variable_id','mz','rt','mz.ms2','rt.ms2','mz_diff','rt_diff','isoform'
    )
)

writexl::write_xlsx(
    list(
        variable_info = object_clean %>% extract_variable_info(),
        expmat = object_clean %>% extract_expression_data() %>% rownames_to_column('variable_id'),
        sample_info = object_clean %>% extract_sample_info()
    ), "RemoveRedundancy/metadata.xlsx"

)

## extract annotation   =====================================
Sys.sleep(3)

message(msg_run('extract annotation'))
future::plan("multicore",workers = 6)
final_anno <- 
furrr::future_map_dfr(1:length(object_neg_list),.f = function(.x){
    tmp.object.neg = object_neg_list[[.x]]
    tmp.neg.anno = tmp.object.neg %>% 
        activate_mass_dataset('variable_info') %>% 
        left_join(neg_variables %>% select(-mz,-rt),by = 'variable_id') %>% 
        filter(tag == 'marked') %>% 
        extract_annotation_table() %>% select(-ms2_files_id)
    
    tmp.object.pos = object_pos_list[[.x]]
    tmp.pos.anno = tmp.object.pos %>% 
        activate_mass_dataset('variable_info') %>% 
        left_join(pos_variables %>% select(-mz,-rt),by = 'variable_id') %>% 
        filter(tag == 'marked') %>% 
        extract_annotation_table() %>% select(-ms2_files_id)

    out = rbind(tmp.neg.anno,tmp.pos.anno)
    return(out)
},  .progress = T)

final_anno_clean <- 
final_anno %>% tibble() %>% 
    arrange(variable_id) %>% 
    filter(
        Adduct == '(M+H)+' | Adduct == "(M-H)-"
    ) %>% 
    group_by(variable_id) %>% 
    filter(Level == min(Level)) %>% 
    ungroup() %>% 
    mutate(
        Compound.name_fix = str_split(Compound.name,";",2,T)[,1]
    ) %>% 
    mutate(
        Compound.name_fix = str_to_lower(Compound.name_fix) %>% str_remove(.," $")
    ) %>% 
    group_by(variable_id,Compound.name_fix) %>% 
    slice_max(Total.score,n = 1) %>% 
    slice_head(n = 1)

writexl::write_xlsx(final_anno_clean,"RemoveRedundancy/Final_annotation.xlsx")
message(msg_yes('Step2 finish'))
