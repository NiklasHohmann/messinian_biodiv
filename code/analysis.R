cat("Running analysis, please wait...\n")

#### load data ####
cat("Loading data\n")
messinian_db <- read.csv(file = "data/messinianDB.csv")

#### Set seet ####
set.seed(1)

#### Replace invalid tax names with NA ####
messinian_db$Species.name = replace(messinian_db$Species.name, messinian_db$Species.name == "sp.", NA)
messinian_db$Genus.name = replace(messinian_db$Genus.name, messinian_db$Genus.name == "indet.", NA)
messinian_db = messinian_db[ !messinian_db$Family == "indet.",  ]

#### Constants ####
timebins <- unique(messinian_db$Age)[c(3, 1, 2)] # sorted from old to young
regions <- unique(messinian_db$region.new)
group.names <- unique(messinian_db$group.name)
group.names.ext = c(group.names, "all groups")
regions.ext = c(regions, "whole basin")
eco_index_names = c("soerensen", "simpson", "nestedness")
timebin_comp = c("T vs. M", "M vs. Z", "T vs. Z")
noOfRep = 10000 # number of repetitions for subsampling

#### Load helper functions ####
cat("Loading helper functions\n")
source("code/helper_functions.R")

#### Species richness through time (whole basin, all groups) ####
cat("Determining species richness\n")
for (group in group.names.ext){
  # extract species names
  Tor = get_from_db(group = group, basin = "whole basin", timeslice = "Tortonian" )
  Mes = get_from_db(group = group, basin = "whole basin", timeslice = "pre-evaporitic Messinian" )
  Zan = get_from_db(group = group, basin = "whole basin", timeslice = "Zanclean" )
  # define subsampling size (80 % of smallest sample)
  subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
  # subsample noOfRep times
  Tor_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
  Mes_sr = rarefyTaxRichness(mySample = Mes, subsampleTo = subsampleTo, noOfRep = noOfRep)
  Zan_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
  # make figure
  file_name = paste0("figs/sr_through_time/sr_through_time_whole_basin_",group,".pdf")
  main = paste0("Species Richness ", group, " whole basin")
  ylim = c(0, max(c(Tor_sr, Mes_sr, Zan_sr)))
  ylab = paste0("Species richness \n subsampled to ",  subsampleTo, " Occurrences")
  pdf(file = file_name)
  boxplot(list( "Tortonian" = Tor_sr,
                "Messinian" = Mes_sr,
                "Zanclean" =  Zan_sr),
          ylim = ylim,
          ylab = ylab,
          main = main)
  dev.off()
}

#### Ecological indices through time (whole basin, all groups) ####
cat("Determining ecological indices\n")
for (group in group.names.ext){
  # extract species names
  Tor = get_from_db(group = group, basin = "whole basin", timeslice = "Tortonian" )
  Mes = get_from_db(group = group, basin = "whole basin", timeslice = "pre-evaporitic Messinian" )
  Zan = get_from_db(group = group, basin = "whole basin", timeslice = "Zanclean" )
  # define subsampling size (80 % of smallest sample)
  subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
  # subsample noOfRep times
  TM = rarefyEcoIndexes(Tor, Mes, subsampleTo, noOfRep)
  MZ = rarefyEcoIndexes(Mes, Zan, subsampleTo, noOfRep)
  TZ = rarefyEcoIndexes(Tor, Zan, subsampleTo, noOfRep)
  # make plots
  for (ind in eco_index_names){
    file_name = paste0("figs/eco_timeslice_comp/",group, "_", ind,"_whole_basin.pdf")
    pdf(file = file_name)
    main = paste0(ind, " for " , group, " whole basin")
    ylab = paste0(ind, "\n subsampled to ",  subsampleTo, " Occurrences")
    boxplot(list("T vs. M" = TM[[ind]],
                 "M vs. Z" = MZ[[ind]],
                 "T vs. Z" = TZ[[ind]]),
            ylim = c(0,1),
            main = main,
            ylab = ylab,
            mar = c(5,5,1,1))
    dev.off()
  }
}

#### Species richness through time for all groups in all regions ####
cat("Determining regional species richness\n")
sr_median = array(data = NA, # mean species richness
                       dim = c(length(group.names.ext), length(timebins), length(regions.ext)),
                       dimnames = list("group" = group.names.ext,
                                       "timebin" = timebins,
                                       "region" = regions.ext))

for (group in group.names.ext){
  for (reg in regions.ext){
    Tor = get_from_db(group = group, basin = reg, timeslice = "Tortonian" )
    Mes = get_from_db(group = group, basin = reg, timeslice = "pre-evaporitic Messinian" )
    Zan = get_from_db(group = group, basin = reg, timeslice = "Zanclean" )
    subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
    Tor_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
    Mes_sr = rarefyTaxRichness(mySample = Mes, subsampleTo = subsampleTo, noOfRep = noOfRep)
    Zan_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
    sr_median[group, "Tortonian", reg] = median(Tor_sr)
    sr_median[group, "pre-evaporitic Messinian", reg] = median(Mes_sr)
    sr_median[group, "Zanclean", reg] = median(Zan_sr)
    
    file_name = paste0("figs/sr_through_time_regional/sr_through_time_", group, "_",  reg ,".pdf")
    ylab = paste0("Species richness \n subsampled to ",  subsampleTo, " Occurrences")
    main = paste0("Species Richness ", group, " ", reg)
    ylim = c(0, max(c(Tor_sr, Mes_sr, Zan_sr)))
    pdf(file = file_name)
    boxplot(list( "Tortonian" = Tor_sr,
                  "Messinian" = Mes_sr,
                  "Zanclean" =  Zan_sr),
            ylim = ylim,
            ylab = ylab,
            main = main)
    dev.off()
    }
}

#### Ecological indices through time per region ( all groups) ####

eco_ind_median = array(data = NA,
                         dim = c(length(group.names.ext), length(timebin_comp), length(regions.ext), length(eco_index_names)),
                         dimnames = list("group" = group.names.ext,
                                         "timebins" = timebin_comp,
                                         "regions" = regions.ext,
                                         "index" = eco_index_names))

cat("Determining regional ecological indices\n")
for (group in c(group.names.ext)){
  for (reg in regions.ext){
    Tor = get_from_db(group = group, basin = reg, timeslice = "Tortonian" )
    Mes = get_from_db(group = group, basin = reg, timeslice = "pre-evaporitic Messinian" )
    Zan = get_from_db(group = group, basin = reg, timeslice = "Zanclean" )
    subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
    TM = rarefyEcoIndexes(Tor, Mes, subsampleTo, noOfRep)
    MZ = rarefyEcoIndexes(Mes, Zan, subsampleTo, noOfRep)
    TZ = rarefyEcoIndexes(Tor, Zan, subsampleTo, noOfRep)
    for (ind in eco_index_names){
      eco_ind_median[group, "T vs. M", reg, ind] = 100 * median(TM[[ind]])
      eco_ind_median[group, "M vs. Z", reg, ind] = 100 * median(MZ[[ind]])
      eco_ind_median[group, "T vs. Z", reg, ind] = 100 * median(TZ[[ind]])
      file_name = paste0("figs/eco_timeslice_comp_regional/",group, "_", reg ,"_", ind, ".pdf")
      pdf(file = file_name)
      main = paste0(ind, " for " , group, " ", reg)
      ylab = paste0(ind, "\n subsampled to ",  subsampleTo, " Occurrences")
      boxplot(list("T vs. M" = TM[[ind]],
                   "M vs. Z" = MZ[[ind]],
                   "T vs. Z" = TZ[[ind]]),
              ylim = c(0,1),
              main = main,
              ylab = ylab,
              mar = c(5,5,1,1))
      dev.off()
    }
  }
}


#### Changes in sr & eco indexes ####
sr_change = array(data = NA,
                  dim = c(length(group.names.ext), 3, length(regions.ext)),
                  dimnames = list("group" = group.names.ext,
                                   "time_slices" = c("TvsM","MvsZ","TvsZ"),
                                  "region" = regions.ext))
for (group in group.names.ext){
  for (reg in regions.ext){
    sr_change[group, "TvsM", reg] = 100 * (1 - sr_median[group, "Tortonian",reg] / sr_median[group,"pre-evaporitic Messinian",reg])
    sr_change[group, "MvsZ", reg] = 100 * (1 - sr_median[group, "pre-evaporitic Messinian",reg] / sr_median[group,"Zanclean",reg])
    sr_change[group, "TvsZ", reg] = 100 * (1 - sr_median[group, "Tortonian",reg] / sr_median[group,"Zanclean",reg])
  }
}

cat("Done. Outputs are in the folder \"figs\" and the variables \"sr_change\" and \"eco_ind_median\".\n")

## assessing sr changes in PERCENT
# first entry: name out of group.names.ext.
# second entry: "TvsM","MvsZ", or "TvsZ": time slice to compare
# third entry: name out of regions.ext: region to examine
# retruns change in median species richness
#sr_change["all groups", "TvsM", "whole basin"] 

## extracting median ecological indices in PERCENT
# first entry: group name, element of group.names.ext
# second entry: time bins to compare, "T vs. M", "M vs. Z" or "T vs. Z"
# third entry: region of interest, element of regions.ext
# fourth entry: ecological index, element of eco_index_names
# returns median ecological index for the specified case.
#eco_ind_median["fish", "T vs. M", "whole basin", "simpson"]


sr_change[,"TvsM", "whole basin"]
sr_change[,"MvsZ", "whole basin"]
eco_ind_median[, "T vs. M", "whole basin", "soerensen"]
eco_ind_median[, "M vs. Z", "whole basin", "soerensen"]
eco_ind_median[, "T vs. Z", "whole basin", "soerensen"]