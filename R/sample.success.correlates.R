rm(list=ls())

library(Mnov.GTseq.data)
library(strataG)
library(dplyr)

data("AS.gtseq")
data("ind.info")
data("Mnov.strata")
data("gtseq.prod.summary")
data("gtseq.val.genos")
data("tissue.archive.data")
data("GTseq.samps.final")
data("id.key")

gtseq.smry.all <- rbind(gtseq.prod.summary, gtseq.val.genos[,1:7])
gtseq.animals <- filter(id.key, LABID %in% gtseq.smry.all$LABID ) %>%
  filter(id.type == "ANIMAL.id") %>% select(c(LABID, ANIMALID = alt.id)) %>%
  left_join(select(ind.info, c("LABID","HAP","SEX")), by = "LABID") %>%
  left_join(select(tissue.archive.data, c("LABID","YR","MO","DA","State","Country","Latitude","Longitude")), by = "LABID") %>%
  left_join(gtseq.smry.all, by = "LABID")

other.sp.smry <- gtseq.smry.all[-which(gtseq.smry.all$LABID %in% gtseq.animals$LABID),] %>%
  left_join(GTseq.samps.final, by = "LABID") %>% 
  select(-c(AnimalID,Duplicate1, Duplicate2, Duplicate3, Herd, HAP, mito.haps, USE, Optimization.sample, Main.run)) 
