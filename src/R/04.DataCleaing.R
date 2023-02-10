##################################
#       Prj: Tidymass
#       Assignment: Data cleaning
#       Author: shawn
#       Date: Feb 9, 2023
#       Location: HENU
##################################
# import packages ---------------------------------------------------------
library(tidymass)
library(tidyverse)
library(MDAtoolkits)
library(IMOtoolkits)


# for negative model ------------------------------------------------------

##> 
# load objects -------------------------------------------------------------------------

setwd("demo/")

load("neg_object")

object.neg <- object

sample_info <- read.delim("injection_order_neg.txt",sep = "\t",header = F)

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
  select(sample_id,group,injection.order,class)

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



# 1.run positive model ---------------------------------------------------------------------

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
  theme(axis.text.x = element_text(size = 3))

ggsave(filename = "01.raw/NEG/missing_value_distribution.png",plot = plt_mv_raw.neg,width = 10,height = 5)
ggsave(filename = "01.raw/NEG/missing_value_distribution.png",plot = plt_mv_raw.neg,width = 10,height = 5)


# 1.2 Outlier detect and remove outliers ----------------------------------

outlier_samples.neg <-
  object.neg %>% 
  `+`(1) %>% 
  log(2) %>% 
  scale() %>% 
  detect_outlier(na_percentage_cutoff = 0.8)

outlier_table.neg <-
  extract_outlier_table(outlier_samples.neg)

out_name.neg <-
  outlier_table.neg %>%
  apply(1, function(x){
    sum(x)
  }) %>%
  `>`(0) %>%
  which()%>% names

if(length(out_name.neg) != 0) {
  
}
