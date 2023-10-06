######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: Add ms2 and annotation
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################

library(tidymass)
library(IMOtoolkits)

# set parameter -----------------------------------------------------------

args<- commandArgs(T)
column <- args[1]

# load objects ------------------------------------------------------------
load("object_neg.rds")
load("object_pos.rds")


# Add MS2 -----------------------------------------------------------------


object_pos.ms2 <-
  object_pos %>%
  mutate_ms2(
    object = .,
    column = column,
    polarity = 'positive',
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "./POS/"
  )
save(object_pos.ms2, file = "object_pos.ms2.rds")

object_neg.ms2 <-
  object_neg %>%
  mutate_ms2(
    object = .,
    column = column,
    polarity = 'negative',
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "./NEG/"
  )
save(object_neg.ms2,file = "object_neg.ms2.rds")



load("~/.HPC_tidymass/MS_db/Agri_plant_KNApSAcK_database.rda")
load("~/.HPC_tidymass/MS_db/RIKEN_PlaSMA_database0.0.1.rda")
load("~/.HPC_tidymass/MS_db/ReSpect_database1.0.rda")
load("~/.HPC_tidymass/MS_db/zma_plantcyc.database.rda")
load("~/.HPC_tidymass/MS_db/mona_database0.0.3.rda")
load("~/.HPC_tidymass/MS_db/massbank_database0.0.3.rda")

# for ms1

object_pos_anno.Agri_knapsack <-
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 7,
    database = knapsack_agri_plant_db,
    candidate.num = 10
  )

object_pos_anno.plantcyc <-
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 7,
    database = zma_plantcyc.database,
    candidate.num = 10
  )


# for ms2

object_pos_anno.mona <-
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 20,
    database = mona_database0.0.3,
    candidate.num = 10
  )

object_pos_anno.massbank <-
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 20,
    database = massbank_database0.0.3,
    candidate.num = 10
  )

object_pos_anno.Plasma <-
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = PlaSMA_database,
    candidate.num = 10
  )

object_pos_anno.respect <-
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = ReSpect_database,
    candidate.num = 10
  )

object_pos_list = 
  c(respect = object_pos_anno.respect,
    plasma = object_pos_anno.Plasma,
    mona = object_pos_anno.mona,
    plantcyc = object_pos_anno.plantcyc,
    knapsack = object_pos_anno.Agri_knapsack)

save(object_pos_list,file = "object_pos_anno.rds")


# neg ---------------------------------------------------------------------

object_neg_anno.Agri_knapsack <-
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 7,
    database = knapsack_agri_plant_db,
    candidate.num = 10
  )

object_neg_anno.plantcyc <-
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 7,
    database = zma_plantcyc.database,
    candidate.num = 10
  )


# for ms2

object_neg_anno.mona <-
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 20,
    database = mona_database0.0.3,
    candidate.num = 10
  )

object_neg_anno.massbank <-
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 20,
    database = massbank_database0.0.3,
    candidate.num = 10
  )

object_neg_anno.Plasma <-
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = PlaSMA_database,
    candidate.num = 10
  )

object_neg_anno.respect <-
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = ReSpect_database,
    candidate.num = 10
  )

object_neg_list = 
  c(respect = object_neg_anno.respect,
   plasma = object_neg_anno.Plasma,
   mona = object_neg_anno.mona,
   plantcyc = object_neg_anno.plantcyc,
   knapsack = object_neg_anno.Agri_knapsack)
save(object_neg_list,file = "object_neg_anno.rds")



