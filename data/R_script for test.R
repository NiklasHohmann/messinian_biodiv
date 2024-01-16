getwd()
setwd("C:/Users/Konstantina/Documents/FISHSALC/Marine biological impact of MSC/review paper/database/FINAL DB/tax validity test against worms")

MessiDB<-read.csv("MessinianDB.csv")
MessiDB$full.name<-paste(MessiDB$Genus.name,MessiDB$Species.name,sep=" ")

library(tidyverse)
library(worrms)

splist<-unique(MessiDB$full.name)

worms_rec <- wm_records_names(name=splist) # returns list of data frames

# bind dataframes into a single dataframe and filter out unaccepted names
worms_rec_MessiDB <- worms_rec %>% bind_rows() %>%
  filter(status=="accepted")
