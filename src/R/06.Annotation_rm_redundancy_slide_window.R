##########################################
#           Prj: HPC-tidymass
#           Assignment: RemoveRedundancy
#           Author: Shawn Wang
#           Date: Wed Jan 03, 2024
##########################################
suppressMessages(library(tidymass))
suppressMessages(library(tidyverse))
suppressMessages(library(furrr))
suppressMessages(library(progressr))
suppressMessages(library(crayon))
msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic
conflicted::conflict_prefer(name = 'map',winner = 'purrr')
handlers(handler_pbcol(
  adjust = 1.0,
  complete = function(s) cli::bg_red(cli::col_black(s)),
  incomplete = function(s) cli::bg_cyan(cli::col_black(s))
))
options(warn = -1)
options(digits = 9)
dir.create("RemoveRedundancy",showWarnings = F,recursive = T)

# data_cleaning -----------------------------------------------------------
message(msg_run('Remove redundancy in positive model....'))
load("object_pos.ms2.rds")
object_pos.ms2
##> check ms1 match ms2
ms2_spectra_info <- 
  object_pos.ms2 %>% 
  extract_ms2_data()
names(ms2_spectra_info) = "QC+Subject"
ms2_spectra_info$`QC+Subject`@variable_id %>% unique() %>% length()
ms2_spectra_info$`QC+Subject`@ms2_spectrum_id %>% unique() %>% length()

##> for mac and linux
conflicted::conflict_prefer('filter','dplyr')
conflicted::conflict_prefer('select','dplyr')

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

mda_check_ms1 = function(x){
  p <- progressor(steps = nrow(x))
  round_n = nrow(x)
  out =  future_map_dfr(.x = 1:nrow(x),.f = function(.x) {
    vid = x[.x,1] |> as.character()
    p(sprintf("%s",paste0('variable_id :',vid,"(",.x,"/",round_n,")")))
    mz1 = x[.x,2] %>% as.numeric()
    rt1 = x[.x,3] %>% as.numeric()
    mz_range = c(mz1 - 20*0.000001,mz1 + 20*0.000001)
    x %>% 
      filter(mz.ms2 > mz_range[1] & mz.ms2 < mz_range[2]) %>% 
      mutate(variable_id_1 = x[.x,1] %>% as.character(),
             delt_mz = (abs(mz1-mz.ms2))*(10^6)/mz1,
             delt_rt = abs(rt1-rt.ms2)) %>% 
      relocate(variable_id_1,.before = variable_id)
  }) %>% 
    filter(variable_id_1 != variable_id) %>% 
    select(-mz_diff,-rt_diff)
  return(out)
}

with_progress({
  check_ms1 <- mda_check_ms1(x = ms1_spec_order_filter_step1)
})

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





candidate_variable_ids <- 
  spec_filter_step2 %>% 
  arrange(desc(rt.ms2)) %>% 
  mutate(rough_rt = round(rt.ms2,0)) %>% 
  group_by(rough_rt) %>% 
  mutate(n = n()) %>% 
  filter(n > 1)


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

rough_rt_index = 
  candidate_variable_ids %>% 
  summarise(n = n()) %>% pull(rough_rt)

vari_list <- map(.x = 1:length(rough_rt_index),.f = function(.x){
  candidate_variable_ids %>% filter(rough_rt == rough_rt_index[.x]) %>% 
    pull(variable_id)
})


prg_run_isf = function(x,rough_rt_index,match_num = 10,n_peak = 20){
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
  x <- x 
  p <- progressr::progressor(along = x);
  round_n = length(x);
  b = furrr::future_map2_dfr(.x = x,.y = c(1:round_n),.f = function(.x,.y){
    rt_g_name = rough_rt_index[.y]
    p(sprintf("%s",paste0('RT-Group:',rt_g_name,"(",.y,"/",round_n,")")))
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
        tmp_out <- tmp_out %>% mutate(RT_group = rt_g_name)
        if(nrow(tmp_out) >= match_num) {
          return(tmp_out )
        } else {
          return(tbl_blank)
        }
      })
    })
    return(tmp.z)
  })
  b <- b %>% drop_na()
  return(b)
}



with_progress({
  res = prg_run_isf(x = vari_list,rough_rt_index = rough_rt_index)
})


res_clean <- 
  res |> group_by(vari_id_1,vari_id_2,RT_group) |> 
  slice_head(n = 1) |> 
  arrange(RT_group) |> 
  as.data.frame()


res_fianl <- 
  res_clean |> 
  select(vari_id_1,vari_id_2,RT_group) |> 
  rowwise() |> 
  mutate(vari_id_1.mz = spec_filter_step2 |> filter(variable_id == vari_id_1) |> pull(mz),
         vari_id_2.mz = spec_filter_step2 |> filter(variable_id == vari_id_2) |> pull(mz),
         mz.diff = vari_id_1.mz - vari_id_2.mz,
         vari_id_1.mean = spec_filter_step2 |> filter(variable_id == vari_id_1) |> pull(mean),
         vari_id_2.mean = spec_filter_step2 |> filter(variable_id == vari_id_2) |> pull(mean),
         main_peak = case_when(
           (abs(mz.diff) < 10) & (vari_id_1.mean > vari_id_2.mean) ~ vari_id_1,
           (abs(mz.diff) < 10) & (vari_id_1.mean < vari_id_2.mean) ~ vari_id_2,
           mz.diff >= 10 ~ vari_id_1,
           mz.diff <= -10 ~ vari_id_2
         ),
         isf_peak= case_when(
           (abs(mz.diff) < 10) & (vari_id_1.mean > vari_id_2.mean) ~ vari_id_2,
           (abs(mz.diff) < 10) & (vari_id_1.mean < vari_id_2.mean) ~ vari_id_1,
           mz.diff >= 10 ~ vari_id_2,
           mz.diff <= -10 ~ vari_id_1
         ))
spec_filter_step3 <- 
  spec_filter_step2 |> 
  anti_join(data.frame(
    variable_id = res_fianl |> pull(isf_peak) |> unique()
  ))|> arrange(mz,rt)


step_loss_num.3  =   nrow(spec_filter_step2) - nrow(spec_filter_step3)
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

prg_run_iso = function(x,rough_mz_index,match_num = 5,n_peak = 10){
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
  x <- x 
  p <- progressr::progressor(along = x);
  round_n = length(x);
  b = furrr::future_map2_dfr(.x = x,.y = c(1:round_n),.f = function(.x,.y){
    mz_g_name = rough_mz_index[.y]
    p(sprintf("%s",paste0('RT-Group:',mz_g_name,"(",.y,"/",round_n,")")))
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
  })
  b <- b %>% drop_na()
  return(b)
}

with_progress({
  res.isoform = prg_run_iso(x = vari_list.iso,rough_mz_index = rough_mz_index,match_num = 10,n_peak = 20)
})

# match_ms2_plt(ms2_spec_order = ms2_spec_order,variable_1 = "M91T260_POS",variable_2 = "M91T223_POS",n_peak = 20)

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
  mutate(tag = "marked") 

writexl::write_xlsx(compound_filted_pos_final,"RemoveRedundancy/Pos_final_rest.xlsx")

# match_ms2_plt(ms2_spec_order = ms2_spec_order,variable_1 = "M80T63_POS",variable_2 = "M81T63_POS",n_peak = 20)


# negative ----------------------------------------------------------------
message(msg_run('Remove redundancy in negative model....'))
load("object_neg.ms2.rds")
object_neg.ms2
##> check ms1 match ms2
ms2_spectra_info <- 
  object_neg.ms2 %>% 
  extract_ms2_data()
names(ms2_spectra_info) = "QC+Subject"
ms2_spectra_info$`QC+Subject`@variable_id %>% unique() %>% length()
ms2_spectra_info$`QC+Subject`@ms2_spectrum_id %>% unique() %>% length()

##> for mac and linux
future::plan("multicore", workers = 7)
conflicted::conflict_prefer('filter','dplyr')
conflicted::conflict_prefer('select','dplyr')

# Step01. multi-match reduandancy -----------------------------------------

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
message(msg_run('Remove redundancy by rough rt range'))
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

##> 2.1 calculate delta mz and delta rt
message(msg_run('Remove redundancy by ms2 spectra....\nmz_tol < 10 ppm & rt_tol < 10 s'))
mda_check_ms1 = function(x){
  p <- progressor(steps = nrow(x))
  round_n = nrow(x)
  out =  future_map_dfr(.x = 1:nrow(x),.f = function(.x) {
    vid = x[.x,1] |> as.character()
    p(sprintf("%s",paste0('variable_id :',vid,"(",.x,"/",round_n,")")))
    mz1 = x[.x,2] %>% as.numeric()
    rt1 = x[.x,3] %>% as.numeric()
    mz_range = c(mz1 - 20*0.000001,mz1 + 20*0.000001)
    x %>% 
      filter(mz.ms2 > mz_range[1] & mz.ms2 < mz_range[2]) %>% 
      mutate(variable_id_1 = x[.x,1] %>% as.character(),
             delt_mz = (abs(mz1-mz.ms2))*(10^6)/mz1,
             delt_rt = abs(rt1-rt.ms2)) %>% 
      relocate(variable_id_1,.before = variable_id)
  }) %>% 
    filter(variable_id_1 != variable_id) %>% 
    select(-mz_diff,-rt_diff)
  return(out)
}

with_progress({
  check_ms1 <- mda_check_ms1(x = ms1_spec_order_filter_step1)
})

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
message(msg_run('Remove ISF.'))
candidate_variable_ids <- 
  spec_filter_step2 %>% 
  arrange(desc(rt.ms2)) %>% 
  mutate(rough_rt = round(rt.ms2,0)) %>% 
  group_by(rough_rt) %>% 
  mutate(n = n()) %>% 
  filter(n > 1)


rough_rt_index = 
  candidate_variable_ids %>% 
  summarise(n = n()) %>% pull(rough_rt)

vari_list <- map(.x = 1:length(rough_rt_index),.f = function(.x){
  candidate_variable_ids %>% filter(rough_rt == rough_rt_index[.x]) %>% 
    pull(variable_id)
})


with_progress({
  res = prg_run_isf(x = vari_list,rough_rt_index = rough_rt_index)
})


res_clean <- 
  res |> group_by(vari_id_1,vari_id_2,RT_group) |> 
  slice_head(n = 1) |> 
  arrange(RT_group) |> 
  as.data.frame()


res_fianl <- 
  res_clean |> 
  select(vari_id_1,vari_id_2,RT_group) |> 
  rowwise() |> 
  mutate(vari_id_1.mz = spec_filter_step2 |> filter(variable_id == vari_id_1) |> pull(mz),
         vari_id_2.mz = spec_filter_step2 |> filter(variable_id == vari_id_2) |> pull(mz),
         mz.diff = vari_id_1.mz - vari_id_2.mz,
         vari_id_1.mean = spec_filter_step2 |> filter(variable_id == vari_id_1) |> pull(mean),
         vari_id_2.mean = spec_filter_step2 |> filter(variable_id == vari_id_2) |> pull(mean),
         main_peak = case_when(
           (abs(mz.diff) < 10) & (vari_id_1.mean > vari_id_2.mean) ~ vari_id_1,
           (abs(mz.diff) < 10) & (vari_id_1.mean < vari_id_2.mean) ~ vari_id_2,
           mz.diff >= 10 ~ vari_id_1,
           mz.diff <= -10 ~ vari_id_2
         ),
         isf_peak= case_when(
           (abs(mz.diff) < 10) & (vari_id_1.mean > vari_id_2.mean) ~ vari_id_2,
           (abs(mz.diff) < 10) & (vari_id_1.mean < vari_id_2.mean) ~ vari_id_1,
           mz.diff >= 10 ~ vari_id_2,
           mz.diff <= -10 ~ vari_id_1
         ))
spec_filter_step3 <- 
  spec_filter_step2 |> 
  anti_join(data.frame(
    variable_id = res_fianl |> pull(isf_peak) |> unique()
  ))|> arrange(mz,rt)


step_loss_num.3  =   nrow(spec_filter_step2) - nrow(spec_filter_step3)
message(
  msg_yes(paste0('Step2. ISF detect  \nA total number of (',step_loss_num.3,"/",nrow(spec_filter_step2),') were removed!'))
)

# isoform identification-----------------------------------------------------------------
message(msg_run('isoform identification'))
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

with_progress({
  res.isoform = prg_run_iso(x = vari_list.iso,rough_mz_index = rough_mz_index,match_num = 10,n_peak = 20)
})


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

# match_ms2_plt(ms2_spec_order = ms2_spec_order,variable_1 = "M229T413_NEG",variable_2 = "M229T444_NEG",n_peak = 20)
message(msg_run('negative model finish, export data.'))
compound_filted_neg_final <- 
spec_filter_step3 %>% 
  left_join(.,iso_label) %>% 
  ungroup() %>% 
  mutate(tag = "marked")

writexl::write_xlsx(compound_filted_neg_final,"RemoveRedundancy/NEG_final_rest.xlsx")

message(msg_yes('Remove redundancy step1 finish!'))


