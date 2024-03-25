##########################################
#           Prj: HPC-tidymass
#           Assignment: RemoveRedundancy
#           Author: Shawn Wang
#           Date: Wed Jan 03, 2024
##########################################
suppressMessages(library(tidymass))
suppressMessages(library(tidyverse))
suppressMessages(library(furrr))
suppressMessages(library(crayon))
msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic
conflicted::conflict_prefer(name = 'map',winner = 'purrr')
##> for mac and linux
conflicted::conflict_prefer('filter','dplyr')
conflicted::conflict_prefer('select','dplyr')
conflicted::conflict_prefer('rename','dplyr')
options(warn = -1)
options(digits = 9)
dir.create("RemoveRedundancy",showWarnings = F,recursive = T)
##> para

step = 1;
windows= 20;
## functions ===================================================================================
##> match ms2 fragment
ms2_fragment_match <- function(ms2_spec_order,variable_1,variable_2,n_peak=20){
  check_index = index %>% 
    filter(variable == variable_1) %>% pull(order) %>% as.numeric()
  check_index2 = index %>% 
    filter(variable == variable_2) %>% pull(order) %>% as.numeric()
  
  check1 = ms2_spec_order[[check_index]] %>% slice_max(order_by = intensity,n = n_peak)
  check2 = ms2_spec_order[[check_index2]] %>% slice_max(order_by = intensity,n = n_peak)
  
  mz_similar <- map_dfr(1:nrow(check1),.f = function(.x){
    frag_mz_1 = check1[.x,1] %>% as.numeric()
    map_dfr(1:nrow(check2),.f = function(.y){
      frag_mz_2 = check2[.y,1] %>% as.numeric()
      mz.diff = abs(frag_mz_1 - frag_mz_2)*(10^6)/frag_mz_1
      tmp.out = data.frame(
        vari_id_1 = variable_1,
        vari_id_2 = variable_2,
        frag_mz_1 = frag_mz_1,
        frag_mz_2 = frag_mz_2,
        frag_mz_1_rank = .x,
        frag_mz_2_rank = .y,
        mz_diff = mz.diff
      )
    }) 
  }) %>% filter(mz_diff < 10) 
  return(mz_similar)
}
##> generate compare group by slide RT windows.
generate_compare_group = function(candidate_group) {
    oplan = future::plan(
    list(
      future::tweak(
        future::multiprocess,
        workers = 8),
      future::tweak(
        future::multicore,
        workers = 10)
    )
  )
  on.exit(future::plan(oplan), add = TRUE)
  temp_out = future_map_dfr(.x = 1:length(candidate_group),.f = function(.x){
    temp_id = candidate_group[[.x]]
    temp_v1_vs_v2 = expand.grid(temp_id,temp_id) %>% distinct() %>% 
      filter(Var1 != Var2 )
  },.progress = T) %>% 
  distinct() %>% 
  rowwise %>% 
  mutate(group_id = sort(c_across(Var1:Var2)) %>% paste0(collapse = '-'))  %>% 
  group_by(group_id) %>% 
  slice_head(n = 1)
  return(temp_out)
}
##> mz and rt diff
mda_check_ms1 = function(x){
  round_n = nrow(x)
  out =  future_map_dfr(.x = 1:nrow(x),.f = function(.x) {
    vid = x[.x,1] |> as.character()
    mz1 = x[.x,2] %>% as.numeric()
    rt1 = x[.x,3] %>% as.numeric()
    mz_range = c(mz1 - 20*0.000001,mz1 + 20*0.000001)
    x %>% 
      filter(mz.ms2 > mz_range[1] & mz.ms2 < mz_range[2]) %>% 
      mutate(variable_id_1 = x[.x,1] %>% as.character(),
             delt_mz = (abs(mz1-mz.ms2))*(10^6)/mz1,
             delt_rt = abs(rt1-rt.ms2)) %>% 
      relocate(variable_id_1,.before = variable_id)
  },.progress = T) %>% 
    filter(variable_id_1 != variable_id) %>% 
    select(-mz_diff,-rt_diff)
  return(out)
}
##> get similar features with has ms2 fragment match at 5/10
get_similar_feature = function(compare_group,n_peak  = 10) {
  temp_var1 = compare_group %>% pull(1) 
  temp_var2 = compare_group %>% pull(2)
      oplan = future::plan(
    list(
      future::tweak(
        future::multiprocess,
        workers = 8),
      future::tweak(
        future::multicore,
        workers = 10)
    )
  )
  on.exit(future::plan(oplan), add = TRUE)
  temp_out = future_map2_dfr(.x = temp_var1,.y = temp_var2,.f = function(.x,.y){
    temp_out_sub = ms2_fragment_match(ms2_spec_order = ms2_spec_order,variable_1 = .x,variable_2 = .y,n_peak = n_peak)
  },.progress = T)
  return(temp_out)
}
##> add_ms_info
get_sim_feature_ms_info = function(similar_features_summary,ms.obj) {
  temp_com_group = similar_features_summary %>% select(1:2,n)
  ##> ms1_info
  vari_info = ms.obj %>% extract_variable_info() %>% select(1:3)
  ##> ms2_info
  vari_info2 = data.frame(
    variable_id = ms2_spectra_obj@variable_id,
    ms2_mz = ms2_spectra_obj@ms2_mz,
    ms2_rt = ms2_spectra_obj@ms2_rt
  )
  tmp_out = temp_com_group %>% 
    inner_join( vari_info %>% setNames(c("vari_id_1","ms1_mz1","ms1_rt1"))) %>% 
    inner_join( vari_info %>% setNames(c("vari_id_2","ms1_mz2","ms1_rt2"))) %>% 
    inner_join( vari_info2 %>% setNames(c("vari_id_1","ms2_mz1","ms2_rt1"))) %>% 
    inner_join( vari_info2 %>% setNames(c("vari_id_2","ms2_mz2","ms2_rt2"))) %>% 
    mutate(
      ms1_mz_diff = abs(ms1_mz1 - ms1_mz2)*(10^6)/ms1_mz2,
      ms2_mz_diff = abs(ms2_mz1 - ms2_mz2)*(10^6)/ms2_mz2,
      ms1_mz_diff_da = abs(ms1_mz1 - ms1_mz2),
      ms2_mz_diff_da = abs(ms2_mz1 - ms2_mz2),
      ms1_rt_diff = abs(ms1_rt1 - ms1_rt2),
      ms2_rt_diff = abs(ms2_rt1 - ms2_rt2)
    )
}
##> check adduct
test_addcut = function(similar_features,add_tbl,mz_tol = 20 ) {
  oplan = future::plan(
    list(
      future::tweak(
        future::multiprocess,
        workers = 8),
      future::tweak(
        future::multicore,
        workers = 10)
    )
  )
  on.exit(future::plan(oplan), add = TRUE)
  out = future_map_dfr(.x = 1:nrow(similar_features),.f = function(.x) {
  temp_add_tbl = add_tbl
  ##> add_cut
  temp_mz_diff = similar_features[.x,12] %>% as.numeric()
  temp_ms2_mz_diff = similar_features[.x,13] %>% as.numeric()
  judge = 
  temp_add_tbl %>% 
    mutate(mz_ms1_diff = abs(abs(mz) - temp_mz_diff)*(10^6)/temp_mz_diff,
           mz_ms2_diff = abs(abs(mz) - temp_ms2_mz_diff)*(10^6)/temp_ms2_mz_diff) %>% 
    filter(mz_ms1_diff < mz_tol | mz_ms2_diff < mz_tol)
  if(nrow(judge) == 0) {
    temp_out = similar_features[.x,] %>% 
      mutate(Adduct = "")
  } else {
    temp_out = similar_features[.x,] %>% 
      mutate(Adduct = judge %>% pull(adduct) %>% paste(.,collapse = ';'))
  }
  return(temp_out)  
},.progress = T) 
  return(out)
}
##> check isf
test_ISF = function(similar_features,ms2_spec_order,n_peak = 20 ) {
  oplan = future::plan(
    list(
      future::tweak(
        future::multiprocess,
        workers = 8),
      future::tweak(
        future::multicore,
        workers = 10)
    )
  )
  on.exit(future::plan(oplan), add = TRUE)

  out = future_map_dfr(.x = 1:nrow(similar_features),function(.x) {
    temp_ms2_mz1 = similar_features[.x,8] %>% as.numeric() %>% setNames(similar_features[.x,1])
    temp_ms2_mz2 = similar_features[.x,10] %>% as.numeric() %>% setNames(similar_features[.x,2])
    temp_ms1_diff = similar_features[.x,12] %>% as.numeric()
    temp_n_fragment = similar_features[.x,3] %>% as.numeric()
    temp_ms1_rt_diff = similar_features[.x,16] %>% as.numeric()
    temp_ms2_rt_diff = similar_features[.x,17] %>% as.numeric()
    temp_ms1_da = similar_features[.x,14] %>% as.numeric()
    temp_ms2_da = similar_features[.x,15] %>% as.numeric()
    if(temp_ms2_mz1 > temp_ms2_mz2) {
      winner = temp_ms2_mz1; loser = temp_ms2_mz2
    } else {
      winner = temp_ms2_mz2; loser = temp_ms2_mz1
    }
    temp_diff = winner - loser
    check_index = index %>% 
      filter(variable == names(winner)) %>% pull(order) %>% as.numeric()
    check1 = ms2_spec_order[[check_index]] %>% slice_max(order_by = intensity,n = n_peak)

    temp_check_isf = 
      check1 %>% 
      mutate(mz_diff_ppm_temp = abs(mz - temp_ms2_mz2)*(10^6)/temp_ms2_mz2,
             frag_diff =  abs(mz - temp_diff)*(10^6)/temp_diff) %>% 
      filter(mz_diff_ppm_temp < 20 | frag_diff < 100)
    
    if(nrow(temp_check_isf) > 0 & (temp_ms1_da > 12 | temp_ms2_da > 12) ) {
      temp_out = similar_features[.x,] %>% 
        mutate(ISF = "type1")
    } else if ((temp_ms1_da > 12 | temp_ms2_da > 12) & (temp_ms1_rt_diff < 2 & temp_ms2_rt_diff < 2)) {
      temp_out = similar_features[.x,] %>% 
        mutate(ISF = "type2")
    } else {
      temp_out = similar_features[.x,] %>% 
      mutate(ISF = "")
    } 
    return(temp_out)
  },.progress = T)
  return(out)
}
##> draw_ms2 
match_ms2_plt = function(ms2_spec_order,variable_1,variable_2,n_peak=20) {
  check_index = index %>% 
    filter(variable == variable_1) %>% pull(order) %>% as.numeric()
  check_index2 = index %>% 
    filter(variable == variable_2) %>% pull(order) %>% as.numeric()
  
  check1 = ms2_spec_order[[check_index]] %>% slice_max(order_by = intensity,n = n_peak)
  check2 = ms2_spec_order[[check_index2]] %>% slice_max(order_by = intensity,n = n_peak)
  check1_df <- 
    check1 %>% 
    mutate(
      variable_id = variable_1,
      x = mz,
      xend = mz,
      y = 0,
      yend = intensity,
      tag = case_when(
        intensity > 10 ~ mz,
        TRUE ~ NA
      )
    )
  check2_df <- 
    check2 %>% 
    mutate(
      variable_id = variable_2 ,
      x = mz,
      xend = mz,
      y = 0,
      yend = -intensity,
      tag = case_when(
        intensity > 10 ~ mz,
        TRUE ~ NA
      )
    )
  
  check_merge = rbind(check1_df,check2_df)
  
  p = ggplot(check_merge,
             mapping = aes(x = x,xend = xend,y = y,yend = yend,color = variable_id)) + 
    geom_segment(size = 1)+
    geom_text(mapping = aes(x = x,y = yend,label = tag),color = 'black')+
    theme_bw() + 
    scale_y_continuous(breaks = c(-100, -50, 0, 50, 100),labels = c('100',"50","0","50","100")) +
    xlab("mz")+
    ylab('intensity') +
    theme(
      axis.ticks = element_line(size = 1),
      panel.border = element_rect(linewidth = 1),
      axis.text = element_text(size = 12)
    )
  return(p)
}
##> isoform
prg_run_iso = function(x,rough_mz_index,match_num = 5,n_peak = 10){
  oplan = future::plan(
    list(
      future::tweak(
        future::multiprocess,
        workers = 4),
      future::tweak(
        future::multicore,
        workers = 6)
    )
  )
  
  on.exit(future::plan(oplan), add = TRUE)
  x <- x 
  round_n = length(x);
  b = furrr::future_map2_dfr(.x = x,.y = c(1:round_n),.f = function(.x,.y){
    mz_g_name = rough_mz_index[.y]
    tbl_blank= data.frame(
      vari_id_1 = NA,
      vari_id_2 = NA,
      frag_mz_1 = NA,
      frag_mz_2 = NA,
      frag_mz_1_rank = NA,
      frag_mz_2_rank = NA,
      mz_diff = NA,
      RT_group = NA
    )
    vari_ids = .x
    n = length(vari_ids)
    tmp.z = map_dfr(1:(n-1),.f = function(.z) {
      tmp.i =  furrr::future_map_dfr((.z + 1):n,.f = function(.i) {
        tmp_out <- ms2_fragment_match(ms2_spec_order = ms2_spec_order,variable_1 = vari_ids[.z],variable_2 = vari_ids[.i],n_peak = n_peak)
        tmp_out <- tmp_out %>% mutate(RT_group = mz_g_name)
        if(nrow(tmp_out) >= match_num) {
          return(tmp_out )
        } else {
          return(tbl_blank)
        }
      })
    })
    return(tmp.z)
  },.progress = T)
  b <- b %>% drop_na()
  return(b)
}


# data_cleaning -----------------------------------------------------------
message(msg_run('Remove redundancy in positive model....'))
load("object_pos.ms2.rds")
obj_qc_mean <- 
object_pos.ms2 %>% extract_expression_data() %>% 
  select(contains('QC')) %>% 
  rownames_to_column("variable_id") %>% 
  group_by(variable_id)%>%
  rowwise() %>%  
  summarise(mean = mean(c_across(contains("QC")))) %>% 
  ungroup()
##> check ms1 match ms2
ms2_spectra_info <- 
  object_pos.ms2 %>% 
  extract_ms2_data()
names(ms2_spectra_info) = "QC+Subject"
ms2_spectra_obj = ms2_spectra_info$`QC+Subject`
ms2_spectra_info$`QC+Subject`@variable_id %>% unique() %>% length()
ms2_spectra_info$`QC+Subject`@ms2_spectrum_id %>% unique() %>% length()

# Step01. multi-match reduandancy -----------------------------------------
message(msg_run('Remove redundancy by rough rt model....'))
##> Get intensity from QC samples
intensity_QC_samples <-
  object_pos.ms2 %>% 
  activate_mass_dataset('expression_data') %>% 
  select(contains('QC')) %>% 
  mutate(mean = rowMeans(.)) %>% 
  extract_expression_data() %>% 
  select(mean) %>% 
  rownames_to_column('variable_id') %>% 
  left_join(object_pos.ms2 %>% extract_variable_info() %>% select(1:3)) %>% 
  arrange(mz,rt)

##> get variable_id
vari_id = ms2_spectra_info$`QC+Subject`@variable_id
index = data.frame(
  variable = vari_id,
  order = 1:length(vari_id)
)

##> extracted features ms1 data
ms2_spec = ms2_spectra_info$`QC+Subject`@ms2_spectra
##> 2.3 reformat ms2 data from massdataset
ms2_spec_order <- map(.x = 1:length(ms2_spec),.f = function(.x){
  ms2_reformed <- 
    ms2_spec[[.x]] %>% as.data.frame() %>% 
    arrange(desc(intensity)) %>% 
    filter(intensity != max(intensity)) %>% 
    mutate(intensity = 100*(intensity/max(intensity)))
  
}) %>% setNames(vari_id)
ms1_spec_order <- map_dfr(.x = 1:length(ms2_spec),.f = function(.x){
  variable_id = vari_id[.x]
  precursor_mz = names(ms2_spec)[.x] %>% str_extract("(?<=mz).*(?=rt)") %>% as.numeric()
  precursor_rt = names(ms2_spec)[.x] %>% str_extract("(?<=rt).*") %>% as.numeric()
  ms2_reformed <- 
    ms2_spec[[.x]] %>% as.data.frame() %>% 
    filter(intensity == max(intensity)) %>% 
    pull(intensity) %>% 
    as.numeric()
  out = data.frame(
    variable_id = variable_id,
    mz.ms2 = precursor_mz,
    rt.ms2 = precursor_rt,
    intensity = ms2_reformed
  )
}) %>% inner_join(intensity_QC_samples) %>% 
  select(variable_id,mz.ms2,rt.ms2,mz,rt,mean) %>% 
  mutate(
    mz_diff = abs(mz - mz.ms2)*(10^6)/mz,
    rt_diff = abs(rt - rt.ms2)
  )


##> rough remove based on rt_diff
ms1_spec_order_filter_step1 <- 
  ms1_spec_order %>% 
  group_by(mz.ms2,rt.ms2) %>% 
  slice_min(rt_diff) %>% 
  slice_head(n = 1)
step_loss_num = nrow(ms1_spec_order) - nrow(ms1_spec_order_filter_step1)
message(
  msg_yes(paste0('Step1. ms1 ms2 spectra multi-match (rough remove, based on rt_diff) \nA total number of (',step_loss_num,"/",nrow(ms1_spec_order),') were removed!'))
)

## based on ms2 info
##> delt mz < 20 ppm
##> delt rt < 30 second
message(msg_run('Remove redundancy by rough ms2 spectra....\nmz_tol < 10 ppm & rt_tol < 10 s'))
##> 2.1 calculate delta mz and delta rt

check_ms1 <- mda_check_ms1(x = ms1_spec_order_filter_step1)



##> 2.2 extract features with delt mz < 10 ppm and rt < 10 second
reduandent_ms1 <- 
  check_ms1 %>% arrange(delt_rt) %>% filter(delt_mz < 10 & delt_rt < 10) %>% 
  as.data.frame() %>%
  rowwise() %>% 
  mutate(group_id = sort(c_across(variable_id_1:variable_id)) %>% paste0(collapse = '-')) %>% 
  group_by(group_id) %>% 
  relocate(group_id,.before = variable_id_1) %>% 
  slice_max(order_by = mean,n = 1)

if(nrow(reduandent_ms1) > 0) {
  ##> 2.4 judge
  ms2_judge <- function(reduandent_ms1,which_row,mz_tol = 10){
    same_ms2_1.index = index %>% filter(variable == reduandent_ms1[which_row,2] %>% as.character()) %>% pull(order)
    same_ms2_2.index = index %>% filter(variable == reduandent_ms1[which_row,3] %>% as.character()) %>% pull(order)
    
    same_ms2_1.ms2 = ms2_spec_order[[same_ms2_1.index]] %>% slice_max(n = 20,order_by = intensity) %>%  pull(mz)
    
    same_ms2_2.ms2 = ms2_spec_order[[same_ms2_2.index]] %>% slice_max(n = 20,order_by = intensity) %>% pull(mz)
    
    tmpx2 <- 
      map_dfr(.x = 1:length(same_ms2_1.ms2),.f = function(.x) {
        map_dfr(1:length(same_ms2_2.ms2),.f = function(.y) {
          data.frame(
            fragment1 = same_ms2_1.ms2[.x],
            fragment2 = same_ms2_2.ms2[.y],
            delt_mz =  abs((same_ms2_1.ms2[.x] - same_ms2_2.ms2[.y])*(10^6)/same_ms2_1.ms2[.x]) 
          )
        })
      }) %>% arrange(delt_mz) %>% filter(delt_mz < mz_tol)
    return(tmpx2)
  }
  
  data_res <- map_dfr(.x = 1:nrow(reduandent_ms1),.f = function(.x) {
    ms2_nrow <- ms2_judge(reduandent_ms1,which_row = .x) %>% nrow()
    out = data.frame(
      group_id = reduandent_ms1[.x,1],
      fragment_align = ms2_nrow
    ) %>% inner_join(reduandent_ms1)
  })
  
  feature_remove_step2 <- 
    data_res %>% select(variable_id)
  
  spec_filter_step2 <- 
    anti_join(
      ms1_spec_order_filter_step1,feature_remove_step2
    ) 
  step_loss_num.2  =   nrow(ms1_spec_order_filter_step1) - nrow(spec_filter_step2)
  message(
    msg_yes(paste0('Step1. ms1 ms2 spectra multi-match (based on ms2 spectra) \nA total number of (',step_loss_num.2,"/",nrow(ms1_spec_order_filter_step1),') were removed!'))
  )
  
} else {
  spec_filter_step2 = ms1_spec_order_filter_step1
  step_loss_num.2  =   nrow(ms1_spec_order_filter_step1) - nrow(spec_filter_step2)
  message(
    msg_yes(paste0('Step1. ms1 ms2 spectra multi-match (based on ms2 spectra) \nA total number of (',step_loss_num.2,"/",nrow(ms1_spec_order_filter_step1),') were removed!'))
  )
}



# In source fragmentation identify ----------------------------------------
message(msg_run('Remove In source fragmentation, ISF'))

# detected adduct ions and isf by RT slide window

spec_filter_step_order <- 
    spec_filter_step2  %>% 
    arrange(rt)

rt_node = c((60+10):(18*60-10))
rt_range_bottom = rt_node - 10
rt_range_top = rt_node + 10


candidate_variable_ids = purrr::map(.x = 1:length(rt_node),.f = function(.x){
  out = spec_filter_step_order  %>% 
  filter(rt > rt_range_bottom[.x] & rt < rt_range_top[.x]) %>% 
  pull(variable_id)
  return(out)
}) %>% unique() 

candidate_single <- candidate_variable_ids %>% keep(.,function(x) length(x) == 1)

candidate_group <- candidate_variable_ids %>% keep(.,function(x) length(x) > 1)

##> generate compare group

compare_group = generate_compare_group(candidate_group = candidate_group)

similar_features <- get_similar_feature(compare_group = compare_group,n_peak = 10)

similar_features_clean <- 
similar_features %>% 
  group_by(vari_id_1,vari_id_2) %>% 
  mutate(n = n()) %>% 
  filter(n >= 5) %>% 
  arrange(desc(n))

similar_features_summary = 
  similar_features_clean %>% 
  filter(frag_mz_1_rank <= 3 & frag_mz_2_rank <= 3) %>% 
  filter(n >= 5) %>% 
  slice_head() %>% 
  arrange(desc(n))

##> add ms info
similar_features_ms_info = get_sim_feature_ms_info(similar_features_summary = similar_features_summary,ms.obj = object_pos.ms2)
##> test adduct
data("rp.pos", envir = environment())
adduct.table <- rp.pos

judge_adduct = test_addcut(similar_features = similar_features_ms_info, add_tbl =adduct.table,mz_tol = 400)
judge_adduct %>% filter(Adduct != "")
##> same_feature 
##> mz_diff < 20 & fragment > 5
similar_features_diff = 
  similar_features_ms_info %>% select(vari_id_1,vari_id_2,n,contains("diff")) %>% 
  inner_join(obj_qc_mean %>% setNames(c("vari_id_1","Intensity_1"))) %>% 
  inner_join(obj_qc_mean %>% setNames(c("vari_id_2","Intensity_2")))

judge_same_feature <- 
  similar_features_diff  %>% 
  filter(
    ms1_mz_diff < 20 | ms2_mz_diff < 20
  )  %>% 
  mutate(
    winner = case_when(
    Intensity_1 > Intensity_2 ~ vari_id_1,
    TRUE ~ vari_id_2
  ),
    loser = case_when(
    Intensity_1 > Intensity_2 ~ vari_id_2,
    TRUE ~ vari_id_1
    ))
same_feature_loser = judge_same_feature %>% pull(loser) %>% unique()
similar_features_step2 = 
  anti_join(
    similar_features_diff,data.frame(vari_id_1 = same_feature_loser)
  ) %>% 
  anti_join(
    data.frame(vari_id_2 = same_feature_loser)
  )

judge_ISF_type2 <- 
  similar_features_step2 %>% 
  filter((ms1_mz_diff_da > 12 | ms2_mz_diff_da > 12) & (ms1_rt_diff < 2 | ms2_rt_diff < 2)) 
variable_info_pos <- 
  object_pos.ms2 %>% extract_variable_info() %>% select(1:3)
ISF_type2_loser <- 
  judge_ISF_type2 %>% 
  left_join(variable_info_pos %>% setNames(c("vari_id_1",'ms1_mz1','ms1_rt1'))) %>% 
  left_join(variable_info_pos %>% setNames(c("vari_id_2",'ms1_mz2','ms1_rt2'))) %>% 
  mutate(
    loser = case_when(
    ms1_mz1 > ms1_mz2 ~ vari_id_2,
    TRUE ~ vari_id_1
    )
  ) %>% pull(loser) %>% unique()

spec_filter_step3 = 
  anti_join(spec_filter_step_order,data.frame(variable_id = c(same_feature_loser,ISF_type2_loser)))

spec_filter_step3 %>% 
  arrange(mz,rt)

step_loss_num.3  =   nrow(spec_filter_step_order) - nrow(spec_filter_step3)
message(
  msg_yes(paste0('Step2. ISF detect  \nA total number of (',step_loss_num.3,"/",nrow(spec_filter_step2),') were removed!'))
)

# isoform identification-----------------------------------------------------------------
message(msg_run('Isoform identification....'))
candidate_variable_ids.iso <- 
  spec_filter_step3 %>% 
  arrange(desc(mz.ms2)) %>% 
  mutate(rough_mz = round(mz.ms2,1)) %>% 
  group_by(rough_mz) %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>% 
  arrange(desc(n))

rough_mz_index = 
  candidate_variable_ids.iso %>% 
  summarise(n = n()) %>% pull(rough_mz)

vari_list.iso <- map(.x = 1:length(rough_mz_index),.f = function(.x){
  candidate_variable_ids.iso %>% filter(rough_mz == rough_mz_index[.x]) %>% 
    pull(variable_id)
})

res.isoform = prg_run_iso(x = vari_list.iso,rough_mz_index = rough_mz_index,match_num = 10,n_peak = 20)


res_clean.iso <- 
  res.isoform |> group_by(vari_id_1,vari_id_2,RT_group) |> 
  slice_head(n = 1) |> 
  arrange(RT_group) |>
  as.data.frame()
iso_label = 
  data.frame(
    variable_id = c(res_clean.iso %>% pull(vari_id_1),res_clean.iso %>% pull(vari_id_2)) %>% unique(),
    isoform = "isoform"
  )
message(msg_run('Positive model finish, export data.'))
compound_filted_pos_final <-  
  spec_filter_step3 %>% 
    left_join(.,iso_label) %>% 
  ungroup() %>% 
  mutate(tag = "marked")  %>% 
  arrange(mz)

writexl::write_xlsx(compound_filted_pos_final,"RemoveRedundancy/POS_final_rest.xlsx")



# negative ----------------------------------------------------------------
message(msg_run('Remove redundancy in negative model....'))
load("object_neg.ms2.rds")
obj_qc_mean <- 
object_neg.ms2 %>% extract_expression_data() %>% 
  select(contains('QC')) %>% 
  rownames_to_column("variable_id") %>% 
  group_by(variable_id)%>%
  rowwise() %>%  
  summarise(mean = mean(c_across(contains("QC")))) %>% 
  ungroup()
##> check ms1 match ms2
ms2_spectra_info <- 
  object_neg.ms2 %>% 
  extract_ms2_data()
names(ms2_spectra_info) = "QC+Subject"
ms2_spectra_obj = ms2_spectra_info$`QC+Subject`
ms2_spectra_info$`QC+Subject`@variable_id %>% unique() %>% length()
ms2_spectra_info$`QC+Subject`@ms2_spectrum_id %>% unique() %>% length()

# Step01. multi-match reduandancy -----------------------------------------
message(msg_run('Remove redundancy by rough rt model....'))
##> Get intensity from QC samples
intensity_QC_samples <-
  object_neg.ms2 %>% 
  activate_mass_dataset('expression_data') %>% 
  select(contains('QC')) %>% 
  mutate(mean = rowMeans(.)) %>% 
  extract_expression_data() %>% 
  select(mean) %>% 
  rownames_to_column('variable_id') %>% 
  left_join(object_neg.ms2 %>% extract_variable_info() %>% select(1:3)) %>% 
  arrange(mz,rt)

##> get variable_id
vari_id = ms2_spectra_info$`QC+Subject`@variable_id
index = data.frame(
  variable = vari_id,
  order = 1:length(vari_id)
)

##> extracted features ms1 data
ms2_spec = ms2_spectra_info$`QC+Subject`@ms2_spectra
##> 2.3 reformat ms2 data from massdataset
ms2_spec_order <- map(.x = 1:length(ms2_spec),.f = function(.x){
  ms2_reformed <- 
    ms2_spec[[.x]] %>% as.data.frame() %>% 
    arrange(desc(intensity)) %>% 
    filter(intensity != max(intensity)) %>% 
    mutate(intensity = 100*(intensity/max(intensity)))
  
}) %>% setNames(vari_id)
ms1_spec_order <- map_dfr(.x = 1:length(ms2_spec),.f = function(.x){
  variable_id = vari_id[.x]
  precursor_mz = names(ms2_spec)[.x] %>% str_extract("(?<=mz).*(?=rt)") %>% as.numeric()
  precursor_rt = names(ms2_spec)[.x] %>% str_extract("(?<=rt).*") %>% as.numeric()
  ms2_reformed <- 
    ms2_spec[[.x]] %>% as.data.frame() %>% 
    filter(intensity == max(intensity)) %>% 
    pull(intensity) %>% 
    as.numeric()
  out = data.frame(
    variable_id = variable_id,
    mz.ms2 = precursor_mz,
    rt.ms2 = precursor_rt,
    intensity = ms2_reformed
  )
}) %>% inner_join(intensity_QC_samples) %>% 
  select(variable_id,mz.ms2,rt.ms2,mz,rt,mean) %>% 
  mutate(
    mz_diff = abs(mz - mz.ms2)*(10^6)/mz,
    rt_diff = abs(rt - rt.ms2)
  )


##> rough remove based on rt_diff
ms1_spec_order_filter_step1 <- 
  ms1_spec_order %>% 
  group_by(mz.ms2,rt.ms2) %>% 
  slice_min(rt_diff) %>% 
  slice_head(n = 1)
step_loss_num = nrow(ms1_spec_order) - nrow(ms1_spec_order_filter_step1)
message(
  msg_yes(paste0('Step1. ms1 ms2 spectra multi-match (rough remove, based on rt_diff) \nA total number of (',step_loss_num,"/",nrow(ms1_spec_order),') were removed!'))
)

## based on ms2 info
##> delt mz < 20 ppm
##> delt rt < 30 second
message(msg_run('Remove redundancy by rough ms2 spectra....\nmz_tol < 10 ppm & rt_tol < 10 s'))
##> 2.1 calculate delta mz and delta rt
check_ms1 <- mda_check_ms1(x = ms1_spec_order_filter_step1)



##> 2.2 extract features with delt mz < 10 ppm and rt < 10 second
reduandent_ms1 <- 
  check_ms1 %>% arrange(delt_rt) %>% filter(delt_mz < 10 & delt_rt < 10) %>% 
  as.data.frame() %>%
  rowwise() %>% 
  mutate(group_id = sort(c_across(variable_id_1:variable_id)) %>% paste0(collapse = '-')) %>% 
  group_by(group_id) %>% 
  relocate(group_id,.before = variable_id_1) %>% 
  slice_max(order_by = mean,n = 1)

if(nrow(reduandent_ms1) > 0) {
  ##> 2.4 judge
  ms2_judge <- function(reduandent_ms1,which_row,mz_tol = 10){
    same_ms2_1.index = index %>% filter(variable == reduandent_ms1[which_row,2] %>% as.character()) %>% pull(order)
    same_ms2_2.index = index %>% filter(variable == reduandent_ms1[which_row,3] %>% as.character()) %>% pull(order)
    
    same_ms2_1.ms2 = ms2_spec_order[[same_ms2_1.index]] %>% slice_max(n = 20,order_by = intensity) %>%  pull(mz)
    
    same_ms2_2.ms2 = ms2_spec_order[[same_ms2_2.index]] %>% slice_max(n = 20,order_by = intensity) %>% pull(mz)
    
    tmpx2 <- 
      map_dfr(.x = 1:length(same_ms2_1.ms2),.f = function(.x) {
        map_dfr(1:length(same_ms2_2.ms2),.f = function(.y) {
          data.frame(
            fragment1 = same_ms2_1.ms2[.x],
            fragment2 = same_ms2_2.ms2[.y],
            delt_mz =  abs((same_ms2_1.ms2[.x] - same_ms2_2.ms2[.y])*(10^6)/same_ms2_1.ms2[.x]) 
          )
        })
      }) %>% arrange(delt_mz) %>% filter(delt_mz < mz_tol)
    return(tmpx2)
  }
  
  data_res <- map_dfr(.x = 1:nrow(reduandent_ms1),.f = function(.x) {
    ms2_nrow <- ms2_judge(reduandent_ms1,which_row = .x) %>% nrow()
    out = data.frame(
      group_id = reduandent_ms1[.x,1],
      fragment_align = ms2_nrow
    ) %>% inner_join(reduandent_ms1)
  })
  
  feature_remove_step2 <- 
    data_res %>% select(variable_id)
  
  spec_filter_step2 <- 
    anti_join(
      ms1_spec_order_filter_step1,feature_remove_step2
    ) 
  step_loss_num.2  =   nrow(ms1_spec_order_filter_step1) - nrow(spec_filter_step2)
  message(
    msg_yes(paste0('Step1. ms1 ms2 spectra multi-match (based on ms2 spectra) \nA total number of (',step_loss_num.2,"/",nrow(ms1_spec_order_filter_step1),') were removed!'))
  )
  
} else {
  spec_filter_step2 = ms1_spec_order_filter_step1
  step_loss_num.2  =   nrow(ms1_spec_order_filter_step1) - nrow(spec_filter_step2)
  message(
    msg_yes(paste0('Step1. ms1 ms2 spectra multi-match (based on ms2 spectra) \nA total number of (',step_loss_num.2,"/",nrow(ms1_spec_order_filter_step1),') were removed!'))
  )
}



# In source fragmentation identify ----------------------------------------
message(msg_run('Remove In source fragmentation, ISF'))

# detected adduct ions and isf by RT slide window

spec_filter_step_order <- 
    spec_filter_step2  %>% 
    arrange(rt)

rt_node = c((60+10):(18*60-10))
rt_range_bottom = rt_node - 10
rt_range_top = rt_node + 10


candidate_variable_ids = purrr::map(.x = 1:length(rt_node),.f = function(.x){
  out = spec_filter_step_order  %>% 
  filter(rt > rt_range_bottom[.x] & rt < rt_range_top[.x]) %>% 
  pull(variable_id)
  return(out)
}) %>% unique() 

candidate_single <- candidate_variable_ids %>% keep(.,function(x) length(x) == 1)

candidate_group <- candidate_variable_ids %>% keep(.,function(x) length(x) > 1)

##> generate compare group

compare_group = generate_compare_group(candidate_group = candidate_group)

similar_features <- get_similar_feature(compare_group = compare_group,n_peak = 10)


similar_features_clean <- 
similar_features %>% 
  group_by(vari_id_1,vari_id_2) %>% 
  mutate(n = n()) %>% 
  filter(n >= 5) %>% 
  arrange(desc(n))

similar_features_summary = 
  similar_features_clean %>% 
  filter(frag_mz_1_rank <= 3 & frag_mz_2_rank <= 3) %>% 
  filter(n >= 5) %>% 
  slice_head() %>% 
  arrange(desc(n))

##> add ms info
similar_features_ms_info = get_sim_feature_ms_info(similar_features_summary = similar_features_summary,ms.obj = object_neg.ms2)
##> test adduct
 data("rp.neg", envir = environment())
adduct.table <- rp.neg

judge_adduct = test_addcut(similar_features = similar_features_ms_info, add_tbl =adduct.table,mz_tol = 400)
judge_adduct %>% filter(Adduct != "")
##> same_feature 
##> mz_diff < 20 & fragment > 5
similar_features_diff = 
  similar_features_ms_info %>% select(vari_id_1,vari_id_2,n,contains("diff")) %>% 
  inner_join(obj_qc_mean %>% setNames(c("vari_id_1","Intensity_1"))) %>% 
  inner_join(obj_qc_mean %>% setNames(c("vari_id_2","Intensity_2")))

judge_same_feature <- 
  similar_features_diff  %>% 
  filter(
    ms1_mz_diff < 20 | ms2_mz_diff < 20
  )  %>% 
  mutate(
    winner = case_when(
    Intensity_1 > Intensity_2 ~ vari_id_1,
    TRUE ~ vari_id_2
  ),
    loser = case_when(
    Intensity_1 > Intensity_2 ~ vari_id_2,
    TRUE ~ vari_id_1
    ))
same_feature_loser = judge_same_feature %>% pull(loser) %>% unique()
similar_features_step2 = 
  anti_join(
    similar_features_diff,data.frame(vari_id_1 = same_feature_loser)
  ) %>% 
  anti_join(
    data.frame(vari_id_2 = same_feature_loser)
  )

judge_ISF_type2 <- 
  similar_features_step2 %>% 
  filter((ms1_mz_diff_da > 12 | ms2_mz_diff_da > 12) & (ms1_rt_diff < 2 | ms2_rt_diff < 2)) 
variable_info_neg <- 
  object_neg.ms2 %>% extract_variable_info() %>% select(1:3)
ISF_type2_loser <- 
  judge_ISF_type2 %>% 
  left_join(variable_info_neg %>% setNames(c("vari_id_1",'ms1_mz1','ms1_rt1'))) %>% 
  left_join(variable_info_neg %>% setNames(c("vari_id_2",'ms1_mz2','ms1_rt2'))) %>% 
  mutate(
    loser = case_when(
    ms1_mz1 > ms1_mz2 ~ vari_id_2,
    TRUE ~ vari_id_1
    )
  ) %>% pull(loser) %>% unique()

spec_filter_step3 = 
  anti_join(spec_filter_step_order,data.frame(variable_id = c(same_feature_loser,ISF_type2_loser)))

spec_filter_step3 %>% 
  arrange(mz,rt)

step_loss_num.3  =   nrow(spec_filter_step_order) - nrow(spec_filter_step3)
message(
  msg_yes(paste0('Step2. ISF detect  \nA total number of (',step_loss_num.3,"/",nrow(spec_filter_step2),') were removed!'))
)

# isoform identification-----------------------------------------------------------------
message(msg_run('Isoform identification....'))
candidate_variable_ids.iso <- 
  spec_filter_step3 %>% 
  arrange(desc(mz.ms2)) %>% 
  mutate(rough_mz = round(mz.ms2,1)) %>% 
  group_by(rough_mz) %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>% 
  arrange(desc(n))

rough_mz_index = 
  candidate_variable_ids.iso %>% 
  summarise(n = n()) %>% pull(rough_mz)

vari_list.iso <- map(.x = 1:length(rough_mz_index),.f = function(.x){
  candidate_variable_ids.iso %>% filter(rough_mz == rough_mz_index[.x]) %>% 
    pull(variable_id)
})
res.isoform = prg_run_iso(x = vari_list.iso,rough_mz_index = rough_mz_index,match_num = 10,n_peak = 20)


res_clean.iso <- 
  res.isoform |> group_by(vari_id_1,vari_id_2,RT_group) |> 
  slice_head(n = 1) |> 
  arrange(RT_group) |>
  as.data.frame()
iso_label = 
  data.frame(
    variable_id = c(res_clean.iso %>% pull(vari_id_1),res_clean.iso %>% pull(vari_id_2)) %>% unique(),
    isoform = "isoform"
  )
message(msg_run('Negative model finish, export data.'))
compound_filted_neg_final <-  
  spec_filter_step3 %>% 
    left_join(.,iso_label) %>% 
  ungroup() %>% 
  mutate(tag = "marked")  %>% 
  arrange(mz)

writexl::write_xlsx(compound_filted_neg_final,"RemoveRedundancy/NEG_final_rest.xlsx")

message(msg_yes('Remove redundancy step1 finish!'))


