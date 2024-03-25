####################################################################################
#           Prj: HPC-tidymass
#           Assignment: semi-annotation
#           Author: Shawn Wang
#           Date: Wed Jan 03, 2024
####################################################################################
suppressMessages(library(tidymass))
suppressMessages(library(tidyverse))
suppressMessages(library(furrr))
suppressMessages(library(crayon))
options(warn = -1)
msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic
dir.create('semi_annotation',showWarnings = F,recursive = T)
load("object_neg.ms2.rds")
load("object_pos.ms2.rds")
variable_final <- readxl::read_xlsx("RemoveRedundancy/metadata_deduplicate_cw_method.xlsx")
## pos   =====================================
message(msg_run('PI and NL calculation in positive model\n'))
vari_ids_pos = variable_final %>% filter(str_detect(variable_id,"POS")) %>% pull(variable_id)
##> extract ms2 data information from massdataset
ms2data_pos = object_pos.ms2 %>% extract_ms2_data() %>% setNames('QC_Subject')

ms2data_pos = ms2data_pos$QC_Subject 

ms2data_pos = ms2data_pos %>% filter(str_detect(variable_id,paste(vari_ids_pos,collapse = "|")))

## PI和NL统计 =====================================
message(msg_run('Extract top 10 fragments of each feature based on fragment ion intensity\n'))
##> intensity 前 10 的碎片
ms2_spc_all = ms2data_pos@ms2_spectra

ms2_PI_pos <-  future_map_dfr(1:length(ms2_spc_all),.f = function(.x) {
    tmp_ms2 <- 
    tmp_fragments <- ms2_spc_all[[.x]] %>% 
        as.data.frame() %>% 
 #       filter(intensity != max(intensity)) %>% 
        slice_max(intensity,n = 10) %>% 
        arrange(desc(intensity)) %>% 
        mutate(mz = round(mz,3),
               variable_id = ms2data_pos@variable_id[.x],
               fragments = paste0('fragment-',str_pad(1:nrow(.),2,'left','0')),
               intensity = (intensity/max(intensity)*100)) %>% 
        select(variable_id,fragments,mz,intensity) 
    return(tmp_ms2)
},.progress = T)
message(msg_run('Calculate NLs by Q1-Qn or Qn-Qn\n'))
ms2_NL_pos <-  future_map_dfr(1:length(ms2_spc_all),.f = function(.x) {
    tmp_fragments <- ms2_spc_all[[.x]] %>% 
        as.data.frame() %>% 
#        filter(intensity != max(intensity)) %>% 
        slice_max(intensity,n = 10) %>% 
        arrange(desc(mz)) %>% 
        pull(mz)
    nls_q1 = future_map(1:length(tmp_fragments),.f = function(.y){
        ms2data_pos@ms2_mz[.x] - tmp_fragments[.y]
    })  %>% unlist()   
    nls_qn = future_map(1:(length(tmp_fragments)-1),.f = function(.y) {
        map_dbl((.y+1):length(tmp_fragments),.f = function(.z) {
            nl = (tmp_fragments[.y] - tmp_fragments[.z]) %>% round(.,3)
        })
    }) %>% unlist()
    nls = c(nls_q1,nls_qn)
    nls = data.frame(
        variabel_id = ms2data_pos@variable_id[.x],
        nls = nls
    )
    return(nls)
},.progress = T) 

ms2_NL_filter_pos <- 
    ms2_NL_pos %>% tibble() %>% 
        filter(nls > 15 & nls < 850) %>% 
        mutate(type = "NL-POS") %>% 
        setNames(c('variable_id','mz','type'))

ms2_PI_filter_pos <- 
    ms2_PI_pos %>% tibble() %>% 
        select(variable_id,mz) %>% 
        filter(mz > 15 & mz < 850) %>% 
        mutate(type = "PI-POS") %>% 
        setNames(c('variable_id','mz','type'))

## neg   =====================================
message(msg_run('PI and NL calculation in negative model\n'))
vari_ids_neg = variable_final %>% filter(str_detect(variable_id,"NEG")) %>% pull(variable_id)
##> extract ms2 data information from massdataset
ms2data_neg = object_neg.ms2 %>% extract_ms2_data() %>% setNames('QC_Subject')

ms2data_neg = ms2data_neg$QC_Subject 

ms2data_neg = ms2data_neg %>% filter(str_detect(variable_id,paste(vari_ids_neg,collapse = "|")))


## PI和NL统计 =====================================
message(msg_run('Extract top 10 fragments of each feature based on fragment ion intensity\n'))
##> intensity 前 10 的碎片
ms2_spc_all = ms2data_neg@ms2_spectra

ms2_PI_neg <-  future_map_dfr(1:length(ms2_spc_all),.f = function(.x) {
    tmp_ms2 <- 
    tmp_fragments <- ms2_spc_all[[.x]] %>% 
        as.data.frame() %>% 
 #       filter(intensity != max(intensity)) %>% 
        slice_max(intensity,n = 10) %>% 
        arrange(desc(intensity)) %>% 
        mutate(mz = round(mz,3),
               variable_id = ms2data_neg@variable_id[.x],
               fragments = paste0('fragment-',str_pad(1:nrow(.),2,'left','0')),
               intensity = (intensity/max(intensity)*100)) %>% 
        select(variable_id,fragments,mz,intensity) 
    return(tmp_ms2)
},.progress = T)
message(msg_run('Calculate NLs by Q1-Qn or Qn-Qn\n'))
ms2_NL_neg <-  future_map_dfr(1:length(ms2_spc_all),.f = function(.x) {
    tmp_fragments <- ms2_spc_all[[.x]] %>% 
        as.data.frame() %>% 
#        filter(intensity != max(intensity)) %>% 
        slice_max(intensity,n = 10) %>% 
        arrange(desc(mz)) %>% 
        pull(mz)
    nls_q1 = future_map(1:length(tmp_fragments),.f = function(.y){
        ms2data_neg@ms2_mz[.x] - tmp_fragments[.y]
    })  %>% unlist()   
    nls_qn = future_map(1:(length(tmp_fragments)-1),.f = function(.y) {
        map_dbl((.y+1):length(tmp_fragments),.f = function(.z) {
            nl = (tmp_fragments[.y] - tmp_fragments[.z]) %>% round(.,3)
        })
    }) %>% unlist()
    nls = c(nls_q1,nls_qn)
    nls = data.frame(
        variabel_id = ms2data_neg@variable_id[.x],
        nls = nls
    )
    return(nls)
},.progress = T) 

ms2_NL_filter_neg <- 
    ms2_NL_neg %>% tibble() %>% 
        filter(nls > 15 & nls < 850) %>% 
        mutate(type = "NL-NEG") %>% 
        setNames(c('variable_id','mz','type'))

ms2_PI_filter_neg <- 
    ms2_PI_neg %>% tibble() %>% 
        select(variable_id,mz) %>% 
        filter(mz > 15 & mz < 850) %>% 
        mutate(type = "PI-NEG") %>% 
        setNames(c('variable_id','mz','type'))
message(msg_run('Merge files\n'))
ms2_fragment_all <- rbind(
    ms2_NL_filter_neg,ms2_NL_filter_pos,ms2_PI_filter_neg,ms2_PI_filter_pos
) 

fragment_summary = 
    ms2_fragment_all %>% 
    group_by(mz,type) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    filter(n >= 3)

message(msg_run('Match with chemical tags based on Chen et.al\n'))
tags <- readxl::read_xlsx("~/.HPC_tidymass/MS_db/Tags_CW.xlsx") %>% 
        mutate(Type_model = case_when(
               str_detect(Fomular,'\\-|\\+') ~ 'PI',
               TRUE ~ "NL" 
        ))
mapped_tags <- 
map_dfr(1:nrow(tags),.f = function(.x){
    Fomular = tags[.x,1] %>% as.character()
    mass_the = tags[.x,2] %>% as.numeric()
    type_model = tags[.x,4] %>% as.character()
    if(type_model == "NL"){
    out = 
    fragment_summary %>% 
        filter(
        (abs(mz-mass_the)*10^6/mz) < 100
        ) %>% 
        mutate(Fomular = Fomular,
               Theoretical_mass = mass_the,
               Type_model = type_model,
               mz_diff = (abs(mz-mass_the)*10^6/mz) 
               ) %>% 
        filter(str_detect(type,"NL"))
    } else if(str_detect(Fomular,"\\-")){
    out = 
    fragment_summary %>% 
        filter(
        (abs(mz-mass_the)*10^6/mz) < 20
        ) %>% 
        mutate(Fomular = Fomular,
               Theoretical_mass = mass_the,
               Type_model = type_model,
               mz_diff = (abs(mz-mass_the)*10^6/mz) 
               ) %>% 
        filter(str_detect(type,"PI-NEG"))
    } else if(str_detect(Fomular,"\\+")) {
    out = 
    fragment_summary %>% 
        filter(
        (abs(mz-mass_the)*10^6/mz) < 20
        ) %>% 
        mutate(Fomular = Fomular,
               Theoretical_mass = mass_the,
               Type_model = type_model,
               mz_diff = (abs(mz-mass_the)*10^6/mz) 
               ) %>% 
        filter(str_detect(type,"PI-POS"))
    }
    return(out)
}) 

semi_anno_final <- 
    inner_join(ms2_fragment_all,mapped_tags) %>%  
    mutate(
        tags = paste0(Type_model,round(Theoretical_mass,0))
    )

semi_anno_final %>% pull(variable_id) %>% unique() %>% length()

semi_anno_final_wider = 
    semi_anno_final %>% 
        select(variable_id,tags) %>% 
        distinct()  %>% 
     group_by(variable_id) %>% 
     summarise(Tag = paste0(tags,collapse = ','))

final_anno <- readxl::read_xlsx("RemoveRedundancy/metadata_deduplicate_cw_method.xlsx",sheet = 4)

message(msg_run('Join semi-annotation and accurate annotation.\n'))
semi_anno_final_tbl <- 
    semi_anno_final_wider %>% 
        left_join(final_anno) %>% 
        select(variable_id,Tag,Lab.ID,Compound.name_fix,Lab.ID,significant)
     
message(msg_run('export labled chemical tags\n'))
features_taged_with_gluco <- 
semi_anno_final %>% filter(str_detect(Fomular,"C6H9O5-|C6H11O6-|C6H10O5|C6H12O6|C6H10O5|C5H8O4|C6H10O4|C5H8O4|C6H10O4")) %>% 
     inner_join(semi_anno_final_tbl) %>% 
     select(variable_id,Tag,Fomular,Lab.ID,Compound.name_fix,significant)  %>% 
     distinct() %>% 
     group_by(variable_id,Lab.ID) %>% 
     mutate(Fomular = paste0(Fomular,collapse = ',')) %>% 
     ungroup() %>% 
     distinct() %>% 
     arrange(significant)

features_taged_with_pheno_amide <- 
    semi_anno_final %>% filter(str_detect(Fomular,"C10H8O3+|C9H8+|C9H7O3+|C9H6O2+|C9H8O2+|C9H10O3+|C9H6O+|C11H10O4+|C9H8O4+|C11H12O3+|C10H12O4+")) %>% 
     inner_join(semi_anno_final_tbl) %>% 
     select(variable_id,Tag,Fomular,Lab.ID,Compound.name_fix,significant)  %>% 
     distinct() %>% 
     group_by(variable_id,Lab.ID) %>% 
     mutate(Fomular = paste0(Fomular,collapse = ',')) %>% 
     ungroup() %>% 
     distinct() %>% 
     arrange(significant)

features_taged_with_oxygen_methyl<- 
    semi_anno_final %>% filter(Fomular == "CH2O") %>% 
     inner_join(semi_anno_final_tbl) %>% 
     select(variable_id,Tag,Fomular,Lab.ID,Compound.name_fix,significant)  %>% 
     distinct() %>% 
     group_by(variable_id,Lab.ID) %>% 
     mutate(Fomular = paste0(Fomular,collapse = ',')) %>% 
     ungroup() %>% 
     distinct() %>% 
     arrange(significant)


features_taged_with_carboxyl<- 
    semi_anno_final %>% filter(Fomular == "CO2") %>% 
     inner_join(semi_anno_final_tbl) %>% 
     select(variable_id,Tag,Fomular,Lab.ID,Compound.name_fix,significant)  %>% 
     distinct() %>% 
     group_by(variable_id,Lab.ID) %>% 
     mutate(Fomular = paste0(Fomular,collapse = ',')) %>% 
     ungroup() %>% 
     distinct() %>% 
     arrange(significant)

features_taged_with_hydroxyl<- 
    semi_anno_final %>% filter(Fomular == "H2O") %>% 
    inner_join(semi_anno_final_tbl) %>% 
        select(variable_id,Tag,Fomular,Lab.ID,Compound.name_fix,significant)  %>% 
        distinct() %>% 
        group_by(variable_id,Lab.ID) %>% 
        mutate(Fomular = paste0(Fomular,collapse = ',')) %>% 
        ungroup() %>% 
        distinct() %>% 
        arrange(significant)


features_taged_with_acetyl<- 
    semi_anno_final %>% filter(Fomular == "C2H2O") %>% 
    inner_join(semi_anno_final_tbl) %>% 
        select(variable_id,Tag,Fomular,Lab.ID,Compound.name_fix,significant)  %>% 
        distinct() %>% 
        group_by(variable_id,Lab.ID) %>% 
        mutate(Fomular = paste0(Fomular,collapse = ',')) %>% 
        ungroup() %>% 
        distinct() %>% 
        arrange(significant)


## export    =====================================

semi_anno_out = list(
    fragment_summary = fragment_summary,
    mapped_tags = mapped_tags,
    semi_anno_final = semi_anno_final,
    semi_anno_final2 = semi_anno_final_tbl,
    acetyl = features_taged_with_acetyl,
    carboxyl = features_taged_with_carboxyl,
    gluco = features_taged_with_gluco,
    hydroxyl = features_taged_with_hydroxyl,
    oxygen_methyl = features_taged_with_oxygen_methyl,
    pheno_amide = features_taged_with_pheno_amide
)
writexl::write_xlsx(semi_anno_out,"semi_annotation/semi_annotation.xlsx")

message(msg_yes('All process finish!\n'))


