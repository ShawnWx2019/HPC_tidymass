######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: Add ms2 and annotation
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################

library(tidymass)
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



load("/home/data/public/01.database/01.Database/02.meta_db/fiehn_hilic_database0.0.3.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/hmdb_database0.0.3.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/kegg_ms1_database0.0.3.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/massbank_database0.0.3.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/mona_database0.0.3.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/snyder_database_hilic0.0.3.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/KNApSAcK_ms1_database.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/plantcyc_ms1_database0.0.2.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/RIKEN_PlaSMA_database0.0.1.rda")
load("/home/data/public/01.database/01.Database/02.meta_db/Natural_products_database_v1.rds")


# for ms1

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos.ms2,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = kegg_ms1_database0.0.3,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = plantcyc_ms1_database0.0.2,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = KNApSAcK_ms1_database0.0.1,
    candidate.num = 2
  )
# for ms2

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = snyder_database_rplc0.0.3,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = mona_database0.0.3,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = massbank_database0.0.3,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = hmdb_database0.0.3,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = fiehn_hilic_database0.0.3,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = RIKEN_PlaSMA_database0.0.1,
    candidate.num = 2
  )

object_pos_anno <-
  annotate_metabolites_mass_dataset(
    object = object_pos_anno,
    polarity = 'positive',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = Natural_products_database_v1,
    candidate.num = 2
  )
save(object_pos_anno,file = "object_pos_anno.rds")


# neg ---------------------------------------------------------------------

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg.ms2,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = kegg_ms1_database0.0.3,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = plantcyc_ms1_database0.0.2,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = KNApSAcK_ms1_database0.0.1,
    candidate.num = 2
  )
# for ms2

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = snyder_database_rplc0.0.3,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = mona_database0.0.3,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = massbank_database0.0.3,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = hmdb_database0.0.3,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = fiehn_hilic_database0.0.3,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = RIKEN_PlaSMA_database0.0.1,
    candidate.num = 2
  )

object_neg_anno <-
  annotate_metabolites_mass_dataset(
    object = object_neg_anno,
    polarity = 'negative',
    ms1.match.ppm = 15,
    column = column,
    threads = 100,
    database = Natural_products_database_v1,
    candidate.num = 2
  )

save(object_neg_anno,file = "object_neg_anno.rds")

