##########################################
#           Prj: HPC-tidymass
#           Assignment: RemoveRedundancy-cw
#           Author: Shawn Wang
#           Date: Wed Jan 03, 2024
##########################################
suppressMessages(library(tidyverse))
suppressMessages(library(tidymass))
suppressMessages(library(furrr))
suppressMessages(library(progressr))
suppressMessages(library(crayon))
options(warn = -1)
msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic
load("object_neg.ms2.rds")
handlers(handler_pbcol(
      adjust = 1.0,
    complete = function(s) cli::bg_red(cli::col_black(s)),
  incomplete = function(s) cli::bg_cyan(cli::col_black(s))
))
## functions   =====================================
##> extract ms2_data from tidymass for next steps
mda_extract_ms2_data = function(object) {
    ms2_data = object %>% extract_ms2_data() %>% setNames('ms2_data')
    ms2_data = ms2_data$ms2_data
    return(ms2_data)
}
##> compare ms2 information of two variables
ms2_compare <- function(ms2_data,variable_1,variable_2, fragments_num = 20,mz_tol = 20) {

    ##> extract ms2 information for each ms2_spectra
    temp1_ms2_data = ms2_data %>% filter(variable_id == variable_1)
    temp2_ms2_data = ms2_data %>% filter(variable_id == variable_2)
    temp1_ms2_id = temp1_ms2_data@ms2_spectrum_id
    temp2_ms2_id = temp2_ms2_data@ms2_spectrum_id
    temp1_mz = temp1_ms2_data@ms2_mz
    temp2_mz = temp2_ms2_data@ms2_mz
    temp1_rt = temp1_ms2_data@ms2_rt
    temp2_rt = temp2_ms2_data@ms2_rt
    temp1_fragments_n = temp1_ms2_data@ms2_spectra %>% .[[1]] %>% as.data.frame() %>% slice_max(intensity,n = fragments_num)
    temp2_fragments_n = temp2_ms2_data@ms2_spectra %>% .[[1]] %>% as.data.frame() %>% slice_max(intensity,n = fragments_num)
    ##> intensity of each Q1 were followed the methods of tidymass. sum of top 5 fragments
    temp1_intensity = temp1_fragments_n  %>% slice_max(intensity,n = 5) %>% pull(intensity) %>% sum()
    temp2_intensity = temp2_fragments_n  %>% slice_max(intensity,n = 5) %>% pull(intensity) %>% sum()
    temp1_fragments_mz = temp1_fragments_n %>% pull(mz)
    temp2_fragments_mz = temp2_fragments_n %>% pull(mz)
    temp_fragments_delta = map(seq_along(temp1_fragments_mz),.f = function(.x){
        map(seq_along(temp2_fragments_mz),.f = function(.y) {
            temp1_ms2_mz = temp1_fragments_mz[.x]
            temp2_ms2_mz = temp2_fragments_mz[.y]
            temp_ms2_mz_diff = abs(temp1_ms2_mz-temp2_ms2_mz)*(10^6)/temp1_ms2_mz
        }) %>% unlist()
    }) %>% unlist()
    temp_fragments_diff = temp_fragments_delta[temp_fragments_delta < mz_tol]
    temp_sign_fragment_num = length(temp_fragments_diff)
    summary = data.frame(
        variable_id_1 = variable_1,
        variable_id_2 = variable_2,
        variable_1_ms2_id = temp1_ms2_id,
        variable_2_ms2_id = temp2_ms2_id,
        variable_1_intensity = temp1_intensity,
        variable_2_intensity = temp2_intensity,
        ms2_mz_diff = abs(temp1_mz - temp2_mz)*(10^6)/temp1_mz,
        ms2_rt_diff = abs(temp1_rt - temp2_rt),
        similar_fragement_num = temp_sign_fragment_num
    )

    return(summary)
}
##> main function for feature deduplicated
mda_cw_deduplicate = function(vari_ids,
                                    features,
                                    ms2_data,
                                    winners = FALSE,
                                    ms1_mz_tol = 20,
                                    ms1_mz_type = 'rel',
                                    rt_tol = 30,
                                    mz_tol = 20,
                                    similar_fragement_num = 5, 
                                    core_num = 80,
                                    fragments_num = 20) {

    if(isFALSE(winners)) {
        ms2_data = ms2_data
        features = features
    } else {
        ms2_data = ms2_data %>% filter(str_detect(variable_id,paste(vari_ids,collapse = "|")))
        features = data.frame(
            variable_id = vari_ids
        ) %>% left_join(features)
    }
    plan("multicore",workers = core_num)
    batch_run = function(vari_ids) {
    p <- progressr::progressor(along = vari_ids);
    round_n = length(vari_ids);
    round1_filter = future_map_dfr(seq_along(vari_ids),.f = function(.x) {
    temp_seed_vari = vari_ids[.x]
    p(sprintf("%s",paste0('variable_id:',temp_seed_vari,"(",.x,"/",round_n,")")))
    temp_seed_mz = features %>% 
        filter(variable_id == temp_seed_vari) %>% 
        pull(mz.ms2)
    if(ms1_mz_type == 'rel'){
        temp_vari_branch = 
        features %>% 
        filter(abs(mz.ms2 - temp_seed_mz)*(10^6)/mz.ms2 < ms1_mz_tol) 
    } else if(ms1_mz_type == 'abs'){
        temp_vari_branch = 
        features %>% 
        filter(abs(mz.ms2 - temp_seed_mz) < ms1_mz_tol) 
    }

    if(nrow(temp_vari_branch) <= 1) {
        out = data.frame(
        variable_id_1 = NA,
        variable_id_2 = NA,
        variable_1_ms2_id = NA,
        variable_2_ms2_id = NA,
        variable_1_intensity = NA,
        variable_2_intensity = NA,
        ms2_mz_diff = NA,
        ms2_rt_diff = NA,
        similar_fragement_num = NA
    )
    } else {
       temp_vari_branch = temp_vari_branch %>% pull(variable_id)
       out = map_dfr(seq_along(temp_vari_branch),.f = function(.y) {
        ms2_compare(ms2_data = ms2_data,
                    variable_1 = temp_seed_vari,
                    variable_2 = temp_vari_branch[.y],
                    fragments_num = fragments_num,
                    mz_tol = mz_tol)
       })
    }
    return(out)  
    }) %>% 
    tibble() %>% 
    drop_na() %>%
    filter(variable_id_1 != variable_id_2)
    }
    with_progress({
        round1_filter = batch_run(vari_ids = vari_ids)
    })
round1_filter_clean <- 
round1_filter %>%  
    filter(ms2_mz_diff < ms1_mz_tol) %>% 
    filter(ms2_rt_diff < rt_tol | similar_fragement_num >= similar_fragement_num ) %>% 
    tibble()  %>% 
    rowwise() %>% 
    mutate(group_id = sort(c_across(variable_id_1:variable_id_2)) %>% paste0(collapse = '-')) %>% 
    relocate(group_id,.before = variable_id_1) %>% 
    group_by(group_id) %>% 
    slice_head(n = 1)
out = 
    round1_filter_clean  %>% 
        select(variable_id_1,variable_id_2,variable_1_intensity,variable_2_intensity) %>% 
        mutate(loser = case_when(
            variable_1_intensity > variable_2_intensity ~ variable_id_2,
            variable_1_intensity <= variable_2_intensity ~ variable_id_1
            ),
            winner = case_when(
            variable_1_intensity > variable_2_intensity ~ variable_id_1,
            variable_1_intensity <= variable_2_intensity ~ variable_id_2)      
        ) %>% 
        ungroup()
return(out)

}
## round1    =====================================
message(msg_run('negative model'))
message(msg_run('Remove redundant features with following cutoff:\nmz tol 20, rt_tol 0.5 fragment_tol 5/20 in 20'))
neg_feature <- readxl::read_xlsx("RemoveRedundancy/Non_redundant_feature.xlsx")
ms2_data = object_neg.ms2 %>% mda_extract_ms2_data()
##> mz tol 20, rt_tol 0.5 fragment_tol 5/20 in 20
vari_ids = neg_feature %>% pull(variable_id)
##> round 1 step1. scan all features
round1_feature_remove = mda_cw_deduplicate(
    vari_ids = vari_ids,
    ms2_data = ms2_data,
    features = neg_feature,
    ms1_mz_tol = 20,
    mz_tol = 20,rt_tol = 30,
    similar_fragement_num = 5,fragments_num = 5,
    core_num = 20
)
round1_loser.tmp =  round1_feature_remove %>% pull(loser) %>% unique() %>% sort()
round1_loser = round1_loser.tmp %>% unique()
vari_ids =  round1_feature_remove %>% pull(winner) %>% unique() %>% sort()
round1_count = nrow(round1_feature_remove)
##> round1 step2. scan winners for loop 
round.num = 1
while(round1_count > 0) {
    round.num = round.num + 1
    message(msg_run(paste0('loop: ',round.num)))
    round1_feature_remove = mda_cw_deduplicate(
        vari_ids = vari_ids,
        ms2_data = ms2_data,
        features = neg_feature,
        ms1_mz_tol = 20,
        mz_tol = 20,rt_tol = 30,
        similar_fragement_num = 5,fragments_num = 5,
        core_num = 20,
        winner = T
    )
    round1_count = nrow(round1_feature_remove)
    if(round1_count == 0) { 
        break ## no more duplicated features were detected.
    } 
    round1_loser.tmp =  round1_feature_remove %>% pull(loser) %>% unique() %>% sort()
    round1_loser = c(round1_loser,round1_loser.tmp) %>% unique() ## merge all features in each round
    vari_ids =  round1_feature_remove %>% pull(winner) %>% unique() %>% sort()
}
message(msg_warning(paste0('After ',round.num,' round, all redundancy features were removed!')))

# library(ggvenn)
# ggvenn(
#     list(winner = vari_ids %>% sort(),
#     loser = round1_loser %>% sort())
# )

##> Remove duplicated features
round1_loser = data.frame(variable_id = round1_loser)
round1_feature_remove_res = 
    anti_join(neg_feature,round1_loser)
nrow(neg_feature)
nrow(round1_feature_remove_res)
## round2 deduplicate    =====================================
message(msg_run('Remove redundant features with following cutoff:\nmz tol 1.5 da, rt_tol 0.3 fragment_tol 2/10 in 10'))
##> extract variables for round2
vari_id_round2 = round1_feature_remove_res %>% pull(variable_id)

round2_feature_remove = mda_cw_deduplicate(
    vari_ids = vari_id_round2,
    features = neg_feature,
    ms2_data = ms2_data,
    rt_tol = 12, mz_tol = 10,
    ms1_mz_tol = 1.5,ms1_mz_type = 'abs',
    similar_fragement_num = 2,fragments_num = 20,
    core_num = 40
)

round2_loser.tmp =  round2_feature_remove %>% pull(loser) %>% unique() %>% sort()
round2_loser = round2_loser.tmp %>% unique()
vari_id_round2 =  round2_feature_remove %>% pull(winner) %>% unique() %>% sort()
round2_count = nrow(round2_feature_remove)

##> round2 step2. scan winners for loop 
round.num = 1
while(round2_count > 0) {
    message(msg_run(paste0('loop: ',round.num)))
    round2_feature_remove = mda_cw_deduplicate(
        vari_ids = vari_id_round2,
        features = neg_feature,
        ms2_data = ms2_data,
        rt_tol = 12, mz_tol = 10,
        ms1_mz_tol = 1.5,ms1_mz_type = 'abs',
        similar_fragement_num = 2,fragments_num = 20,
        core_num = 40,winners = T
    )

    round2_loser.tmp =  round2_feature_remove %>% pull(loser) %>% unique() %>% sort()
    round2_loser = c(round2_loser,round2_loser.tmp) %>% unique() ## merge all features in each round
    vari_id_round2 =  round2_feature_remove %>% pull(winner) %>% unique() %>% sort()
    round2_count = nrow(round2_feature_remove)
}
message(msg_warning(paste0('After ',round.num,' round, all redundancy features were removed!')))

##> Remove duplicated features
round2_loser = data.frame(variable_id = round2_loser)
round2_feature_remove_res = 
    anti_join(round1_feature_remove_res,round2_loser)

feature_final_negative = round2_feature_remove_res
nrow(round2_feature_remove_res)



## positive    =====================================
message(msg_run('positive model'))
message(msg_run('Remove redundant features with following cutoff:\nmz tol 20, rt_tol 0.5 fragment_tol 5/20 in 20'))
load("object_pos.ms2.rds")

## round1    =====================================
pos_feature <- readxl::read_xlsx("RemoveRedundancy/Non_redundant_feature.xlsx",sheet = 2)
ms2_data = object_pos.ms2 %>% mda_extract_ms2_data()
##> mz tol 20, rt_tol 0.5 fragment_tol 5/20 in 20
vari_ids = pos_feature %>% pull(variable_id)
##> round 1 step1. scan all features
round1_feature_remove = mda_cw_deduplicate(
    vari_ids = vari_ids,
    ms2_data = ms2_data,
    features = pos_feature,
    ms1_mz_tol = 20,
    mz_tol = 20,rt_tol = 30,
    similar_fragement_num = 5,fragments_num = 5,
    core_num = 20
)
round1_loser.tmp =  round1_feature_remove %>% pull(loser) %>% unique() %>% sort()
round1_loser = round1_loser.tmp %>% unique()
vari_ids =  round1_feature_remove %>% pull(winner) %>% unique() %>% sort()
round1_count = nrow(round1_feature_remove)
##> round1 step2. scan winners for loop 
round.num = 1
while(round1_count > 0) {
    round.num = round.num + 1
    message(msg_run(paste0('loop: ',round.num)))
    round1_feature_remove = mda_cw_deduplicate(
        vari_ids = vari_ids,
        ms2_data = ms2_data,
        features = pos_feature,
        ms1_mz_tol = 20,
        mz_tol = 20,rt_tol = 30,
        similar_fragement_num = 5,fragments_num = 5,
        core_num = 20,
        winner = T
    )
    round1_count = nrow(round1_feature_remove)
    if(round1_count == 0) { 
        break ## no more duplicated features were detected.
    } 
    round1_loser.tmp =  round1_feature_remove %>% pull(loser) %>% unique() %>% sort()
    round1_loser = c(round1_loser,round1_loser.tmp) %>% unique() ## merge all features in each round
    vari_ids =  round1_feature_remove %>% pull(winner) %>% unique() %>% sort()
}
message(msg_warning(paste0('After ',round.num,' round, all redundancy features were removed!')))

# library(ggvenn)
# ggvenn(
#     list(winner = vari_ids %>% sort(),
#     loser = round1_loser %>% sort())
# )

##> Remove duplicated features
round1_loser = data.frame(variable_id = round1_loser)
round1_feature_remove_res = 
    anti_join(pos_feature,round1_loser)

## round2 deduplicate    =====================================
message(msg_run('Remove redundant features with following cutoff:\nmz tol 1.5 da, rt_tol 0.3 fragment_tol 2/10 in 10'))
##> extract variables for round2
vari_id_round2 = round1_feature_remove_res %>% pull(variable_id)

round2_feature_remove = mda_cw_deduplicate(
    vari_ids = vari_id_round2,
    features = pos_feature,
    ms2_data = ms2_data,
    rt_tol = 12, mz_tol = 10,
    ms1_mz_tol = 1.5,ms1_mz_type = 'abs',
    similar_fragement_num = 2,fragments_num = 20,
    core_num = 40
)

round2_loser.tmp =  round2_feature_remove %>% pull(loser) %>% unique() %>% sort()
round2_loser = round2_loser.tmp %>% unique()
vari_id_round2 =  round2_feature_remove %>% pull(winner) %>% unique() %>% sort()
round2_count = nrow(round2_feature_remove)
round.num = 1
##> round2 step2. scan winners for loop 
while(round2_count > 0) {
    round.num = round.num + 1
    message(msg_run(paste0('loop: ',round.num)))
    round2_feature_remove = mda_cw_deduplicate(
        vari_ids = vari_id_round2,
        features = pos_feature,
        ms2_data = ms2_data,
        rt_tol = 12, mz_tol = 10,
        ms1_mz_tol = 1.5,ms1_mz_type = 'abs',
        similar_fragement_num = 2,fragments_num = 20,
        core_num = 40,winners = T
    )

    round2_loser.tmp =  round2_feature_remove %>% pull(loser) %>% unique() %>% sort()
    round2_loser = c(round2_loser,round2_loser.tmp) %>% unique() ## merge all features in each round
    vari_id_round2 =  round2_feature_remove %>% pull(winner) %>% unique() %>% sort()
    round2_count = nrow(round2_feature_remove)
}
message(msg_warning(paste0('After ',round.num,' round, all redundancy features were removed!')))
##> Remove duplicated features
round2_loser = data.frame(variable_id = round2_loser)
round2_feature_remove_res = 
    anti_join(round1_feature_remove_res,round2_loser)


feature_final_positive = round2_feature_remove_res
feature_anno <- readxl::read_xlsx("RemoveRedundancy/Final_annotation.xlsx")
feature_anno <- 
    feature_anno %>% select(variable_id,Lab.ID,Compound.name_fix,Compound.name,KEGG.ID,Total.score)
feature_final_all = rbind(feature_final_negative,feature_final_positive) 
feature_anno_all <- dplyr::inner_join(x = feature_final_all %>% select(variable_id),y = feature_anno, by = 'variable_id') %>% 
                    distinct() %>% 
    mutate(significant = case_when(
        str_detect(Lab.ID,'ReSpect|PlaSMA') ~ 1,
        str_detect(Lab.ID,'MONA') ~ 2,
        str_detect(Lab.ID,'Zma') ~ 3,
        TRUE ~ 4
    ))

feature_anno_represent <- 
    feature_anno_all %>% 
    group_by(variable_id) %>% 
    slice_min(significant) %>% 
    slice_max(Total.score)

feature_anno_represent %>% 
    pull(variable_id) %>% 
    unique() %>% length()

expmat <-readxl::read_xlsx("RemoveRedundancy/metadata.xlsx",sheet = 2)
expmat_clean <- 
    expmat %>% select(-starts_with("QC")) %>% 
    inner_join(
        feature_final_all %>% select(variable_id)
    )
message(msg_run("export metadata"))
writexl::write_xlsx(list(
    feature_infomation = feature_final_all,
    expmat = expmat_clean,
    annotaiton_all = feature_anno_all,
    feature_anno_represent = feature_anno_represent
),"RemoveRedundancy/metadata_deduplicate_cw_method.xlsx")
message(msg_yes("remove redundancy by Chen method finish!"))
