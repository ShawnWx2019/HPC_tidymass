##################################
#       Prj: Tidymass
#       Assignment: Data cleaning
#       Author: shawn
#       Date: Feb 9, 2023
#       Location: HENU
##################################
# import packages ---------------------------------------------------------
library(tidymass)
suppressMessages(if (!require('devtools')) BiocManager::install('devtools'))
suppressMessages(if (!require('tidyverse')) BiocManager::install('tidyverse'))
suppressMessages(if (!require('patchwork')) BiocManager::install('patchwork'))
suppressMessages(if (!require('MDAtoolkits')) install_github(repo = "ShawnWx2019/MDAtoolkits",ref = 'master'))
suppressMessages(if (!require('IMOtoolkits')) install_github(repo = "ShawnWx2019/IMOtoolkits"))
# library(tidymass)
# for negative model ------------------------------------------------------
args <- commandArgs(T)

sample_info_path <- args[2]
heterogeneous <- args[1]
# load objects -------------------------------------------------------------------------


load("object.neg")

object.neg <- object

sample_info <- read.delim("injection_order_neg.txt",sep = "\t",header = F)
sample_info2 <- read.delim(sample_info_path,sep = "\t",header = T)
sample_info <-
  sample_info %>% 
  mutate(
    sample_id = str_extract(
      string = V1,
      pattern = "(?<= )\\w*\\d(?=.raw)"
    ),
    class = case_when(
      str_detect(sample_id,"QC") ~ "QC",
      str_detect(sample_id,"S") ~ "Subject"
    )
  ) %>% 
  drop_na() %>% 
  mutate(injection.order = c(1:nrow(.))) %>% 
  select(sample_id,injection.order,class) %>% 
  left_join(sample_info2,by = "sample_id") %>% 
  mutate(
    group = case_when(
      class == "QC" ~ "QC",
      TRUE ~ group
    )
  )

# function ----------------------------------------------------------------

##> theme

theme1 =   theme(
  panel.border = element_rect(size = 1.5),
  axis.title = element_text(size = 14,color = 'black'),
  axis.text = element_text(size = 12,color = 'black')
)

##> add sample info
add_sample_info = function(object, sample_info = sample_info) {
  object <- 
    object %>% 
    activate_mass_dataset('sample_info') %>% 
    dplyr::select('sample_id') %>% 
    left_join(sample_info,"sample_id")
  return(object)
}

##> remove batch effect
batch_detect = function(object) {
  plt_batch =
    object %>%
    extract_expression_data() %>%
    dplyr::select(contains("QC")) %>%
    pivot_longer(contains("QC"),values_to = "value",names_to = "QC_samples") %>%
    drop_na() %>%
    mutate(
      value = log2(value),
      QC_samples = paste0(
        "QC",str_pad(
          string = gsub("QC","",QC_samples),width = 3,side = 'left',pad = "0"
        )
      )
    ) %>%
    arrange(QC_samples) %>%
    ggplot(data = .,mapping = aes(x = QC_samples,y = value,color = QC_samples))+
    geom_boxplot()+
    ylab("log2(Raw peak area)")+
    xlab("")+
    theme_bw()+
    theme1+
    theme(
      axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),
      legend.position = "none"
    )
  return(plt_batch)
}

## peak distribution
peak_distribution <- function(object) {
  plt_peak_dis =
    object %>%
    `+`(1) %>%
    log(10) %>%
    show_mz_rt_plot() +
    scale_size_continuous(range = c(0.01, 2))+
    theme1
  return(plt_peak_dis)
}



# 1.run negative model ---------------------------------------------------------------------

object.neg = add_sample_info(
  object = object.neg,
  sample_info = sample_info
)  

##> save raw object
##> 

dir.create("01.raw/NEG",showWarnings = F,recursive = T)
save(object.neg,file = "01.raw/NEG/object.neg")


# 1.1 Overview of raw data  ---------------------------------------------

##> Detected batch effect.
batch_detect_plt.neg <- 
  batch_detect(object = object.neg)

ggsave(filename = "01.raw/NEG/raw_peak_boxplot.neg.png",plot = batch_detect_plt.neg,width = 10,height = 5)
ggsave(filename = "01.raw/NEG/raw_peak_boxplot.neg.pdf",plot = batch_detect_plt.neg,width = 10,height = 5)

##> peak distribution
##> 

peak_distribution_plt.neg <-
  peak_distribution(object = object.neg)
ggsave(filename = "01.raw/NEG/peak_distribution_plt.neg.png",plot = peak_distribution_plt.neg,width = 10,height = 5)
ggsave(filename = "01.raw/NEG/peak_distribution_plt.neg.pdf",plot = peak_distribution_plt.neg,width = 10,height = 5)

##> missing value
plt_mv_raw.neg<-
  show_sample_missing_values(object = object.neg, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "01.raw/NEG/missing_value_distribution.png",plot = plt_mv_raw.neg,width = 10,height = 5)
ggsave(filename = "01.raw/NEG/missing_value_distribution.pdf",plot = plt_mv_raw.neg,width = 10,height = 5)


# 1.2 Outlier detect and remove outliers ----------------------------------

qc_id = object.neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id <- 
  object.neg %>% 
  activate_mass_dataset(
    what = "sample_info"
  ) %>% 
  dplyr::filter(class == "Subject") %>% 
  pull(sample_id)

object.neg <-
  object.neg %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>% 
  mutate_variable_na_freq(according_to_samples = subject_id) 
##> na frequence is less than 0.2 in qc samples.
if(heterogeneous=="no"){
  object.neg.mv <-
    object.neg %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(na_freq < 0.2 & na_freq.1 < 0.5)
} else {
  object.neg.mv <-
    object.neg %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(na_freq < 0.2)
}

dir.create("02.remove_noise/NEG",showWarnings = F,recursive = T)
save(object.neg.mv,file = "02.remove_noise/NEG/object.neg.mv")

plt_mv_remove_noise.neg<-
  show_sample_missing_values(object = object.neg.mv, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 8,angle = 90)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "02.remove_noise//NEG/plt_mv_remove_noise.neg.png",plot = plt_mv_remove_noise.neg,width = 10,height = 5)
ggsave(filename = "02.remove_noise//NEG/plt_mv_remove_noise.neg.pdf",plot = plt_mv_remove_noise.neg,width = 10,height = 5)

if(heterogeneous=="no") {
  outlier_samples.neg <-
    object.neg.mv %>% 
    `+`(1) %>% 
    log(2) %>% 
    scale() %>% 
    detect_outlier(na_percentage_cutoff = 0.8)
  
  outlier_table.neg <-
    extract_outlier_table(outlier_samples.neg)

  out_name.neg <-
    outlier_table.neg %>%
    filter(according_to_na == TRUE) %>% 
    rownames()
  ##> remove outlier based on 
  if(length(out_name.neg) != 0) {
    object.neg.mv <- 
      object.neg.mv %>% 
      activate_mass_dataset('expression_data') %>% 
      select(-all_of(out_name.neg))
  }
}


plt_mv_remove_outlier.neg<-
  show_sample_missing_values(object = object.neg.mv, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 8,angle = 90)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "02.remove_noise//NEG/plt_mv_remove_outlier.neg.png",plot = plt_mv_remove_outlier.neg,width = 10,height = 5)
ggsave(filename = "02.remove_noise//NEG/plt_mv_remove_outlier.neg.pdf",plot = plt_mv_remove_outlier.neg,width = 10,height = 5)


# impute mv ---------------------------------------------------------------
object.neg.impute =
  object.neg.mv %>%
  impute_mv(method = 'knn')
dir.create("03.impute_mv/NEG",showWarnings = F,recursive = T)
save(object.neg.impute,file = "03.impute_mv/NEG/object.neg.impute")


# normalization -----------------------------------------------------------
dir.create("04.normalization/NEG",recursive = T,showWarnings = F)
object_norm_neg <-
  normalize_data(object = object.neg.impute,method = "svr",threads = 80)

object_inte_neg <-
  integrate_data(object_norm_neg,'qc_mean')
save(object_inte_neg,file = "04.normalization/NEG/object_inte_neg")
##> pca_plot before normalization
pca_before =
  IMO_plt_pca(
    object.neg.impute,tag = "group",scale = T,interactive = F,showImage = T
  )
##> pca plot after normalization 
pca_after <- 
  IMO_plt_pca(
    object_inte_neg,tag = "group",scale = T,interactive = F,showImage = T
  )

plt_pca_BandA_norm = (pca_before$plot+ggtitle("before"))+
  theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 2)))+
  (pca_after$plot+ggtitle('after'))+
  theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 2)))+
  plot_layout(guides = "collect")+
  plot_annotation(
    theme = theme(legend.position = "bottom"),tag_levels = "A"
  )

ggsave(plot = plt_pca_BandA_norm,filename = "04.normalization/NEG/plt_pca_BandA_norm.png",width = 11,height = 7)
ggsave(plot = plt_pca_BandA_norm,filename = "04.normalization/NEG/plt_pca_BandA_norm.pdf",width = 11,height = 7)

##> check rsd.
plt_rsd_BandA <-
  IMO_plt_rsd(
    obj_old = object.neg.impute,
    obj_new = object_inte_neg,QC_tag = "QC",interactive = F,showImage = T,x_loc = c(130,130),y_loc = c(15,100)
    )
ggsave(filename = "04.normalization/NEG/plt_rsd_BandA.png",width = 8,height = 7)
ggsave(filename = "04.normalization/NEG/plt_rsd_BandA.pdf",width = 8,height = 7)

##> remove features with large rsd
neg_rsd = plt_rsd_BandA$rsd_tbl
object_neg <- object_inte_neg %>% 
  activate_mass_dataset('variable_info') %>% 
  left_join(neg_rsd,by = c('variable_id' = 'ID') ) %>% 
  filter(norm.rsd <= 30)
save(object_neg,file = "object_neg.rds")

# Positive model ----------------------------------------------------------

load("object.pos")

object.pos <- object

sample_info <- read.delim("injection_order_pos.txt",sep = "\t",header = F)

sample_info <-
  sample_info %>% 
  mutate(
    sample_id = str_extract(
      string = V1,
      pattern = "(?<= )\\w*\\d(?=.raw)"
    ),
    group = case_when(
      str_detect(sample_id,"QC") ~ "QC",
      str_detect(sample_id,"S") ~ "Subject"
    ),
    class = group
  ) %>% 
  drop_na() %>% 
  mutate(injection.order = c(1:nrow(.))) %>% 
  select(sample_id,injection.order,class) %>% 
  left_join(sample_info2,by = "sample_id") %>% 
  mutate(
    group = case_when(
      class == "QC" ~ "QC",
      TRUE ~ group
    )
  )



# 2.run positive model ---------------------------------------------------------------------

object.pos = add_sample_info(
  object = object.pos,
  sample_info = sample_info
)  

##> save raw object
##> 

dir.create("01.raw/POS",showWarnings = F,recursive = T)
save(object.pos,file = "01.raw/POS/object.pos")


# 1.1 Overview of raw data  ---------------------------------------------

##> Detected batch effect.
batch_detect_plt.pos <- 
  batch_detect(object = object.pos)

ggsave(filename = "01.raw/POS/raw_peak_boxplot.pos.png",plot = batch_detect_plt.pos,width = 10,height = 5)
ggsave(filename = "01.raw/POS/raw_peak_boxplot.pos.pdf",plot = batch_detect_plt.pos,width = 10,height = 5)

##> peak distribution
##> 

peak_distribution_plt.pos <-
  peak_distribution(object = object.pos)
ggsave(filename = "01.raw/POS/peak_distribution_plt.pos.png",plot = peak_distribution_plt.pos,width = 10,height = 5)
ggsave(filename = "01.raw/POS/peak_distribution_plt.pos.pdf",plot = peak_distribution_plt.pos,width = 10,height = 5)

##> missing value
plt_mv_raw.pos<-
  show_sample_missing_values(object = object.pos, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "01.raw/POS/missing_value_distribution.png",plot = plt_mv_raw.pos,width = 10,height = 5)
ggsave(filename = "01.raw/POS/missing_value_distribution.pdf",plot = plt_mv_raw.pos,width = 10,height = 5)


# 1.2 Outlier detect and remove outliers ----------------------------------


qc_id = object.pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id <- 
  object.pos %>% 
  activate_mass_dataset(
    what = "sample_info"
  ) %>% 
  dplyr::filter(class == "Subject") %>% 
  pull(sample_id)

object.pos <-
  object.pos %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>% 
  mutate_variable_na_freq(according_to_samples = subject_id)

##> na frequence is less than 0.2 in qc samples.
##> 

if(heterogeneous=="no"){
  object.pos.mv <-
    object.pos %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(na_freq < 0.2 & na_freq.1 < 0.5)
} else {
  object.pos.mv <-
    object.pos %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(na_freq < 0.2)
}

dir.create("02.remove_noise/POS",showWarnings = F,recursive = T)
save(object.pos.mv,file = "02.remove_noise/POS/object.pos.mv")

plt_mv_remove_noise.pos<-
  show_sample_missing_values(object = object.pos.mv, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 8,angle = 90)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "02.remove_noise/POS/plt_mv_remove_noise.pos.png",plot = plt_mv_remove_noise.pos,width = 10,height = 5)
ggsave(filename = "02.remove_noise/POS/plt_mv_remove_noise.pos.pdf",plot = plt_mv_remove_noise.pos,width = 10,height = 5)

##> remove outlier

if(heterogeneous=="no") {
  outlier_samples.pos <-
    object.pos.mv %>% 
    `+`(1) %>% 
    log(2) %>% 
    scale() %>% 
    detect_outlier(na_percentage_cutoff = 0.5)
  
  outlier_table.pos <-
    extract_outlier_table(outlier_samples.pos)
  
  out_name.pos <-
    outlier_table.pos %>%
    filter(according_to_na == TRUE) %>% 
    rownames()
  
  ##> remove outlier based on 
  if(length(out_name.pos) != 0) {
    object.pos.mv <- 
      object.pos.mv %>% 
      activate_mass_dataset('expression_data') %>% 
      select(-all_of(out_name.pos))
  }
}


plt_mv_remove_outlier.pos<-
  show_sample_missing_values(object = object.pos.mv, percentage = TRUE,color_by = 'group',order_by = 'injection.order')+theme1+
  geom_hline(yintercept = 80,color = "red",linetype="dashed")+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 8,angle = 90)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()+
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "02.remove_noise/POS/plt_mv_remove_outlier.pos.png",plot = plt_mv_remove_outlier.pos,width = 10,height = 5)
ggsave(filename = "02.remove_noise/POS/plt_mv_remove_outlier.pos.pdf",plot = plt_mv_remove_outlier.pos,width = 10,height = 5)


# impute mv ---------------------------------------------------------------
object.pos.impute =
  object.pos.mv %>%
  impute_mv(method = 'knn')
dir.create("03.impute_mv/POS",showWarnings = F,recursive = T)
save(object.pos.impute,file = "03.impute_mv/POS/object.pos.impute")


# normalization -----------------------------------------------------------
dir.create("04.normalization/POS",recursive = T,showWarnings = F)
object_norm_pos <-
  normalize_data(object = object.pos.impute,method = "svr",threads = 80)

object_inte_pos <-
  integrate_data(object_norm_pos,'qc_mean')
save(object_inte_pos,file = "04.normalization/POS/object_inte_pos")
##> pca_plot before normalization
pca_before =
  IMO_plt_pca(
    object.pos.impute,tag = "group",scale = T,interactive = F,showImage = T
  )
##> pca plot after normalization 
pca_after <- 
  IMO_plt_pca(
    object_inte_pos,tag = "group",scale = T,interactive = F,showImage = T
  )

plt_pca_BandA_norm = (pca_before$plot+ggtitle("before"))+
  theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 2)))+
  (pca_after$plot+ggtitle('after'))+
  theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 2)))+
  plot_layout(guides = "collect")+
  plot_annotation(
    theme = theme(legend.position = "bottom"),tag_levels = "A"
  )

ggsave(plot = plt_pca_BandA_norm,filename = "04.normalization/POS/plt_pca_BandA_norm.png",width = 11,height = 7)
ggsave(plot = plt_pca_BandA_norm,filename = "04.normalization/POS/plt_pca_BandA_norm.pdf",width = 11,height = 7)

##> check rsd.
plt_rsd_BandA <-
  IMO_plt_rsd(
    obj_old = object.pos.impute,
    obj_new = object_inte_pos,QC_tag = "QC",interactive = F,showImage = T,x_loc = c(130,130),y_loc = c(15,100)
  )
ggsave(filename = "04.normalization/POS/plt_rsd_BandA.png",width = 8,height = 7)
ggsave(filename = "04.normalization/POS/plt_rsd_BandA.pdf",width = 8,height = 7)

##> remove features with large rsd
pos_rsd = plt_rsd_BandA$rsd_tbl
object_pos <- object_inte_pos %>% 
  activate_mass_dataset('variable_info') %>% 
  left_join(pos_rsd,by = c('variable_id' = 'ID')) %>% 
  filter(norm.rsd <= 30)

save(object_pos,file = "object_pos.rds")
