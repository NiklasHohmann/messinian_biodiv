#### load data ####
messinian_db <- read.csv(file = "data/messinianDB.csv")

#### Replace invalid tax names with NA ####
messinian_db$Species.name = replace(messinian_db$Species.name, messinian_db$Species.name == "sp.", NA)
messinian_db$Genus.name = replace(messinian_db$Genus.name, messinian_db$Genus.name == "indet.", NA)
messinian_db = messinian_db[ !messinian_db$Family == "indet.",  ]

#### Constants ####
timebins <- unique(messinian_db$Age)[c(3, 1, 2)] # sorted from old to young
regions <- unique(messinian_db$region.new)
group.names <- unique(messinian_db$group.name)
noOfRep = 10000

#### Load helper functions ####
source("code/helper_functions.R")

#### Species richness through time (whole basin, all groups) ####
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


# #### eco indexes geographic comparison through time ####
# ## get sample isze
# samplesize=c()
# for (timebin in timebins){
#   for (basin in regions){
#     samplesize=c(samplesize,length(get_from_db(taxLevel='species',group='all groups',basin=basin,timeslice=timebin)))
#   }
# }
# 
# subsampleTo=ceiling(0.8*min(samplesize))
# n=5000
# for (timeslice in timebins){
#   Eas=get_from_db(taxLevel='species',
#                 group="all groups",
#                 basin="Eastern Mediterranean",
#                 timeslice=timeslice)
#   PoP=get_from_db(taxLevel='species',
#                 group="all groups",
#                 basin="Po Plain-Northern Adriatic",
#                 timeslice=timeslice)
#   Wes=get_from_db(taxLevel='species',
#                 group="all groups",
#                 basin="Western Mediterranean",
#                 timeslice=timeslice)
#   EasVsPoP=rarefyEcoIndexes(Eas,PoP,subsampleTo,noOfRep = n)
#   PoPVsWes=rarefyEcoIndexes(PoP,Wes,subsampleTo,noOfRep = n)
#   WesVsEas=rarefyEcoIndexes(Wes,Eas,subsampleTo,noOfRep = n)
#   pdf(paste("figs/", "Soerensen " ,timeslice ,", all groups, Species level, comparison between regions.pdf"))
#   boxplot(list("East vs. Po Plain"=EasVsPoP$soerensen,
#                "Po Plain vs. West"=PoPVsWes$soerensen,
#                "West vs. East"=WesVsEas$soerensen),
#           main=paste("Soerensen ", timeslice ,", all groups, Species level, comparison between regions"),
#           ylab=paste("Subsampled to", subsampleTo, "occurrences"))
#   dev.off()
#   
#   pdf(paste("figs/", "Simpson " ,timeslice ,", all groups, Species level, comparison between regions.pdf"))
#   boxplot(list("East vs. Po Plain"=EasVsPoP$simpson,
#                "Po Plain vs. West"=PoPVsWes$simpson,
#                "West vs. East"=WesVsEas$simpson),
#           main=paste("Simpson ", timeslice ,", all groups, Species level, comparison between regions"),
#           ylab=paste("Subsampled to", subsampleTo, "occurrences"))
#   dev.off()
#   
#   pdf(paste("figs/", "Nestedness " ,timeslice ,", all groups, Species level, comparison between regions.pdf"))
#   boxplot(list("East vs. Po Plain"=EasVsPoP$nestedness,
#                "Po Plain vs. West"=PoPVsWes$nestedness,
#                "West vs. East"=WesVsEas$nestedness),
#           main=paste("Nestedness ", timeslice ,", all groups, Species level, comparison between regions"),
#           ylab=paste("Subsampled to", subsampleTo, "occurrences"))
#   dev.off()
#   print(timeslice)
# }
# 
# 
# #### species richenss
# regions=c("Eastern Mediterranean","Po Plain-Northern Adriatic","Western Mediterranean")
# for (timeslice in timebins){
#   Eas=get_from_db(taxLevel='species',
#                 group="all groups",
#                 basin="Eastern Mediterranean",
#                 timeslice=timeslice)
#   PoP=get_from_db(taxLevel='species',
#                 group="all groups",
#                 basin="Po Plain-Northern Adriatic",
#                 timeslice=timeslice)
#   Wes=get_from_db(taxLevel='species',
#                 group="all groups",
#                 basin="Western Mediterranean",
#                 timeslice=timeslice)
#   subsampleTo=900 #floor(max(1,0.8*min(c(length(Eas),length(PoP),length(Wes)))))
#   print(subsampleTo)
#   speciesRichnessList=list("Eastern Mediterranean"=rarefyTaxRichness(Eas,
#                                                                      subsampleTo = subsampleTo,
#                                                                      noOfRep = 50000),
#                            "Po Plain-Northern Adriatic"=rarefyTaxRichness(PoP,
#                                                                           subsampleTo = subsampleTo,
#                                                                           noOfRep = 50000),
#                            "Western Mediterranean"=rarefyTaxRichness(Wes,
#                                                                      subsampleTo = subsampleTo,
#                                                                      noOfRep = 50000))
#   ymax=800 #max(sapply(speciesRichnessList,function(x) max(x)))
#   pdf(paste("figs/", "Species_richness " ,timeslice ,", all groups, Species level.pdf"))
#   boxplot(speciesRichnessList,
#           main=paste("Species Richness ", timeslice ,", all groups, Species level"),
#           ylab=paste("Subsampled to", subsampleTo, "occurrences"),
#           ylim=c(0,1.1*ymax))
#   dev.off()
#   
# }
# 
# 
# rarefyTaxRichness(a,subsampleTo = 10,noOfRep = 1000)
# 
# group.names=sort(unique(messiDB$group.name),na.last = TRUE)
# 
# #### fill in missing genus names with family names
# has.genus.name=!is.na(messiDB$Genus.name)
# messiDB$Genus.name[!has.genus.name]=paste(messiDB$Family[!has.genus.name]," indet.", sep='')
# 
# all(!is.na(messiDB$Genus.name))
# n=10000
# 
# messiDB$Species.name[!is.na(messiDB$Species.name)]= paste(messiDB$Genus.name[!is.na(messiDB$Species.name)], messiDB$Species.name[!is.na(messiDB$Species.name)], sep=' ')
# 
# 
# dataTo=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$Age==timebins[1]]
# dataPre=messiDB$Species.name[!is.na(messiDB$Species.name)  & messiDB$Age==timebins[2]]
# dataZa=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$Age==timebins[3]]
# 
# minN=min(length(dataTo),length(dataPre),length(dataZa))
# sampleSize=max(1,floor(0.8*minN))
# print(sampleSize)
# 
# TOME=myrarefaction(dataTo,dataPre,sampleSize,noofrep = n)
# MEZA=myrarefaction(dataPre,dataZa,sampleSize,noofrep = n)
# 
# ymax=max(c(TOME$sr1,TOME$sr2,MEZA$sr2))
# pdf(paste("figs/", "Species_richness (Whole Basin), all groups, Species level.pdf"))
# boxplot(list("Tortonian"=TOME$sr1,
#              "pre-evaporitic Messinian"=TOME$sr2,
#              "Zanclean"=MEZA$sr2 ),
#         main=paste("Species Richness (Whole Basin) , all groups, Species level"),
#         ylab=paste("Subsampled to", sampleSize, "occurrences"),
#         ylim=c(0,1.1*ymax))
# dev.off()
# 
# pdf(paste("figs/", "Soerensen (Whole Basin)  all groups Species level.pdf"))
# boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$soerensen,
#              "pre-evo Messinian\nvs.\nZanclean"=MEZA$soerensen),
#         main="Soerensen (Whole Basin)   all groups  Species level",
#         ylab=paste("Subsampled to", sampleSize, "occurrences"),
#         ylim=c(0,1))
# dev.off()
# 
# pdf(paste("figs/", "Simpson (Whole Basin)   all groups   Species level.pdf"))
# boxplot(list("Tortonian\nvs.\npre-evo Messinian"=TOME$simpson,
#              "pre-evo Messinian\nvs.\n Zanclean"=MEZA$simpson),
#         main="Simpson (Whole Basin)   all groups   Species level",
#         ylab=paste("Subsampled to", sampleSize, "occurrences"),
#         ylim=c(0,1))
# dev.off()
# 
# pdf(paste("figs/", "Nestedness (Whole Basin)  all groups  species level .pdf"))
# boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$nestedness,
#              "pre-evo Messinian\nvs.\n Zanclean"=MEZA$nestedness),
#         main="Nestedness (Whole Basin)   all groups  Species level",
#         ylab=paste("Subsampled to", sampleSize, "occurrences"),
#         ylim=c(0,1))
# dev.off()
# 
# #### genus level ####
# dataTo=messiDB$Genus.name[!is.na(messiDB$Genus.name) & messiDB$Age==timebins[1]]
# dataPre=messiDB$Genus.name[!is.na(messiDB$Genus.name)  & messiDB$Age==timebins[2]]
# dataZa=messiDB$Genus.name[!is.na(messiDB$Genus.name) & messiDB$Age==timebins[3]]
# 
# minN=min(length(dataTo),length(dataPre),length(dataZa))
# sampleSize=max(1,floor(0.8*minN))
# print(sampleSize)
# 
# TOME=myrarefaction(dataTo,dataPre,sampleSize,noofrep = n)
# MEZA=myrarefaction(dataPre,dataZa,sampleSize,noofrep = n)
# 
# ymax=max(c(TOME$sr1,TOME$sr2,MEZA$sr2))
# pdf(paste("figs/", "Genus_richness (Whole Basin), all groups, Genus level.pdf"))
# boxplot(list("Tortonian"=TOME$sr1,
#              "pre-evaporitic Messinian"=TOME$sr2,
#              "Zanclean"=MEZA$sr2 ),
#         main=paste("Genus Richness (Whole Basin) , all groups, Genus level"),
#         ylab=paste("Subsampled to", sampleSize, "occurrences"),
#         ylim=c(0,1.1*ymax))
# dev.off()
# 
# 
# #####
# 
# for ( group.name in group.names){
#   # whole basin #
#     dataTo=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$group.name == group.name & messiDB$Age==timebins[1]]
#     dataPre=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$group.name == group.name & messiDB$Age==timebins[2]]
#     dataZa=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$group.name == group.name & messiDB$Age==timebins[3]]
#     minN=min(length(dataTo),length(dataPre),length(dataZa))
#     sampleSize=max(1,floor(0.8*minN))
#     print(sampleSize)
#     
#     TOME=myrarefaction(dataTo,dataPre,sampleSize,noofrep = n)
#     MEZA=myrarefaction(dataPre,dataZa,sampleSize,noofrep = n)
#     
#     ymax=max(c(TOME$sr1,TOME$sr2,MEZA$sr2))
#     pdf(paste("figs/", "Species_richness (Whole Basin)", group.name," Species level.pdf"))
#     boxplot(list("Tortonian"=TOME$sr1,
#                  "pre-evaporitic Messinian"=TOME$sr2,
#                  "Zanclean"=MEZA$sr2 ),
#             main=paste("Species Richness (Whole Basin) ", group.name," Species level"),
#             ylab=paste("Subsampled to", sampleSize, "occurrences"),
#             ylim=c(0,1.1*ymax))
#     dev.off()
#     
#     pdf(paste("figs/", "Soerensen (Whole Basin)  ", group.name," Species level.pdf"))
#     boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$soerensen,
#                  "pre-evo Messinian\nvs.\nZanclean"=MEZA$soerensen),
#             main="Soerensen (Whole Basin)  ", group.name," Species level",
#             ylab=paste("Subsampled to", sampleSize, "occurrences"),
#             ylim=c(0,1))
#     dev.off()
#     
#     pdf(paste("figs/", "Simpson (Whole Basin)  ", group.name,"  Species level.pdf"))
#     boxplot(list("Tortonian\nvs.\npre-evo Messinian"=TOME$simpson,
#                  "pre-evo Messinian\nvs.\n Zanclean"=MEZA$simpson),
#             main="Simpson (Whole Basin)  ", group.name,"  Species level",
#             ylab=paste("Subsampled to", sampleSize, "occurrences"),
#             ylim=c(0,1))
#     dev.off()
#     
#     pdf(paste("figs/", "Nestedness (Whole Basin) ", group.name,"species level .pdf"))
#     boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$nestedness,
#                  "pre-evo Messinian\nvs.\n Zanclean"=MEZA$nestedness),
#             main="Nestedness (Whole Basin)  ", group.name," Species level",
#             ylab=paste("Subsampled to", sampleSize, "occurrences"),
#             ylim=c(0,1))
#     dev.off()
#     
#     print(group.name)
# }
# 
# tax.level=c("Family", "Genus.name","Species.name")
# rd=list()
# names=c("fish")
# for(i in 1:length(names)){
#   occ.name=paste(names[i],"_occ",sep='')
#   loc.name=paste(names[i],"_loc",sep='')
#   rd[[names[i]]]$occ=read.csv(paste(occ.name,".csv",sep=''),stringsAsFactors=FALSE)
#   rd[[names[i]]]$loc=read.csv(paste(loc.name,".csv",sep=''),stringsAsFactors=FALSE)
#   rd$fish$occ$Geographic.distribution=sapply(rd$fish$occ$Locality, function(x) rd$fish$loc$Region[x==rd$fish$loc$Locality])
# }
# 
# #regional, sorted by times slices
# rare_res1=list()
# for (timebin in timebins){
#   rare_res1[[timebin]]=list()
#   for (tax in tax.level){
#     dataE=rd$fish$occ[tax][rd$fish$occ$Age==timebin & rd$fish$occ$Geographic.distribution=="Eastern Mediterranean",]
#     dataW=rd$fish$occ[tax][rd$fish$occ$Age==timebin & rd$fish$occ$Geographic.distribution=="Western Mediterranean",]
#     dataP=rd$fish$occ[tax][rd$fish$occ$Age==timebin & rd$fish$occ$Geographic.distribution=="Po Plain-Northern Adriatic",]
#     sampleSize=max(1,floor(0.8*min(length(dataE),length(dataW),length(dataP))))
# 
#     EW=myrarefaction(dataE,dataW,sampleSize,noofrep=1000)
#     WP=myrarefaction(dataW,dataP,sampleSize,noofrep=1000)
#     PE=myrarefaction(dataP,dataE,sampleSize,noofrep=1000)
#     #rare_res1[timebin][tax]=list()
#     rare_res1[[timebin]][[tax]]$EW=EW
#     rare_res1[[timebin]][[tax]]$WP=WP
#     rare_res1[[timebin]][[tax]]$PE=PE
#     rare_res1[[timebin]][[tax]]$nosamples=sampleSize
#     jpeg(paste("figs/", "00fish species richness",tax, timebin, "comparison regions.jpeg"))
#     boxplot(list(Eastern=rare_res1[[timebin]][[tax]]$EW$sr1,
#                  Western=rare_res1[[timebin]][[tax]]$WP$sr1,
#                  Poplain=rare_res1[[timebin]][[tax]]$PE$sr1),
#             main=paste("species richness",tax, timebin))
#     dev.off()
#     
#     jpeg(paste("figs/", "00fish soerensen",tax, timebin, "comparison regions.jpeg"))
#     boxplot(list(EvsW=rare_res1[[timebin]][[tax]]$EW$soerensen,
#                  WvsP=rare_res1[[timebin]][[tax]]$WP$soerensen,
#                  PvsE=rare_res1[[timebin]][[tax]]$PE$soerensen),
#             main=paste("Soerensen",tax, timebin))
#     dev.off()
#     
#     jpeg(paste("figs/", "00fish simpson",tax, timebin, "comparison regions.jpeg"))
#     boxplot(list(EvsW=rare_res1[[timebin]][[tax]]$EW$simpson,
#                  WvsP=rare_res1[[timebin]][[tax]]$WP$simpson,
#                  PvsE=rare_res1[[timebin]][[tax]]$PE$simpson),
#             main=paste("Simpson",tax, timebin))
#     dev.off()
#     
#     jpeg(paste("figs/", "00fish nestedness",tax, timebin, "comparison regions.jpeg"))
#     boxplot(list(EvsW=rare_res1[[timebin]][[tax]]$EW$nestedness,
#                  WvsP=rare_res1[[timebin]][[tax]]$WP$nestedness,
#                  PvsE=rare_res1[[timebin]][[tax]]$PE$nestedness),
#             main=paste("Nestedness",tax, timebin))
#     dev.off()
#   }
# }
# 
# 
# 
# 
# 
# 
# # throug time, sorted by regions
# rare_res2=list()
# for (region in regions){
#   rare_res2[[region]]=list()
#   for (tax in tax.level){
#     dataZ=rd$fish$occ[tax][rd$fish$occ$Age=="Zanclean" & rd$fish$occ$Geographic.distribution==region,]
#     dataM=rd$fish$occ[tax][rd$fish$occ$Age=="pre-evaporitic Messinian" & rd$fish$occ$Geographic.distribution==region,]
#     dataT=rd$fish$occ[tax][rd$fish$occ$Age=="Tortonian" & rd$fish$occ$Geographic.distribution==region,]
# 
#     sampleSize=max(1,floor(0.8*min(length(dataZ),length(dataM),length(dataT))))
#     
#     ZM=myrarefaction(dataZ,dataM,sampleSize,noofrep=1000)
#     MT=myrarefaction(dataM,dataT,sampleSize,noofrep=1000)
#     rare_res2[[region]][[tax]]$ZM=ZM
#     rare_res2[[region]][[tax]]$MT=MT
#     
#     jpeg(paste("figs/", "00fish species richness",tax, region, "comparison time.jpeg"))
#     boxplot(list(Tortonian=rare_res2[[region]][[tax]]$MT$sr2,
#                  peMessinian=rare_res2[[region]][[tax]]$MT$sr1,
#                  Zanclean=rare_res2[[region]][[tax]]$ZM$sr1 ),
#             main=paste("species richness",tax, region))
#     dev.off()
#     
#     jpeg(paste("figs/", "00fish soerensen",tax, region, "comparison time.jpeg"))
#     boxplot(list(TortonianvsppeMessinian=rare_res2[[region]][[tax]]$MT$soerensen,
#                  peMessinianvsZanclean=rare_res2[[region]][[tax]]$ZM$soerensen ),
#             main=paste("soerensen",tax, region))
#     dev.off()
#     
#     jpeg(paste("figs/", "00fish simpson",tax, region, "comparison time.jpeg"))
#     boxplot(list(TortonianvsppeMessinian=rare_res2[[region]][[tax]]$MT$simpson,
#                  peMessinianvsZanclean=rare_res2[[region]][[tax]]$ZM$simpson ),
#             main=paste("simpson",tax, region))
#     dev.off()
#     
#     jpeg(paste("figs/", "00fish nestedness",tax, region, "comparison time.jpeg"))
#     boxplot(list(TortonianvsppeMessinian=rare_res2[[region]][[tax]]$MT$nestedness,
#                  peMessinianvsZanclean=rare_res2[[region]][[tax]]$ZM$nestedness ),
#             main=paste("nestedness",tax, region))
#     dev.off()
#     
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
