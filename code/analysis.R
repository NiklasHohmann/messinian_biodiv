cat("Running analysis, please wait...\n")

#### load data ####
cat("Loading data\n")
messinian_db <- read.csv(file = "data/messinianDB.csv")

#### Replace invalid tax names with NA ####
messinian_db$Species.name = replace(messinian_db$Species.name, messinian_db$Species.name == "sp.", NA)
messinian_db$Genus.name = replace(messinian_db$Genus.name, messinian_db$Genus.name == "indet.", NA)
messinian_db = messinian_db[ !messinian_db$Family == "indet.",  ]

#### Constants ####
timebins <- unique(messinian_db$Age)[c(3, 1, 2)] # sorted from old to young
regions <- unique(messinian_db$region.new)
group.names <- unique(messinian_db$group.name)
noOfRep = 10000 # number of repetitions for subsampling

#### Load helper functions ####
cat("Loading helper functions\n")
source("code/helper_functions.R")

#### Species richness through time (whole basin, all groups) ####
cat("Determining species richness\n")
for (group in c(group.names, "all groups")){
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
for (group in c(group.names, "all groups")){
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
  for (ind in c("soerensen", "simpson", "nestedness")){
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
for (group in c(group.names, "all groups")){
  for (reg in regions){
    Tor = get_from_db(group = group, basin = reg, timeslice = "Tortonian" )
    Mes = get_from_db(group = group, basin = reg, timeslice = "pre-evaporitic Messinian" )
    Zan = get_from_db(group = group, basin = reg, timeslice = "Zanclean" )
    subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
    Tor_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
    Mes_sr = rarefyTaxRichness(mySample = Mes, subsampleTo = subsampleTo, noOfRep = noOfRep)
    Zan_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
    file_name = paste0("figs/sr_through_time_regional/sr_through_time_", group, "_",  reg ,".pdf")
    ylab = "Species Richness"
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
cat("Determining regional ecological indices\n")
for (group in c(group.names, "all groups")){
  for (reg in regions){
    Tor = get_from_db(group = group, basin = "whole basin", timeslice = "Tortonian" )
    Mes = get_from_db(group = group, basin = "whole basin", timeslice = "pre-evaporitic Messinian" )
    Zan = get_from_db(group = group, basin = "whole basin", timeslice = "Zanclean" )
    subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
    TM = rarefyEcoIndexes(Tor, Mes, subsampleTo, noOfRep)
    MZ = rarefyEcoIndexes(Mes, Zan, subsampleTo, noOfRep)
    TZ = rarefyEcoIndexes(Tor, Zan, subsampleTo, noOfRep)
    for (ind in c("soerensen", "simpson", "nestedness")){
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
cat("Done. Outputs are in the folder \"figs\".\n")