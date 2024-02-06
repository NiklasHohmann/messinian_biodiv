# load DB first

timebins=c("Tortonian", "pre-evaporitic Messinian","Zanclean")
regions=c("Eastern Mediterranean","Po Plain-Northern Adriatic","Western Mediterranean")
group.names= c("benthic_foraminifera","bivalves", "bryozoans" ,  "corals",  "dinoflagellates","echinoids","fish", "gastropods" ,"marine_mammals","nanoplankton", "ostracods" ,  "planktic_foraminifera", "scaphopod_chitons_cephalopods", "sharks") 

#### Aux Functions
{
MSCSample=function(taxLevel,group,basin,timeslice){
  #### groups 
  stopifnot(group %in% c('all groups',group.names))
  if (group=='all groups'){
    groupIndex=rep(TRUE,length(messiDB$group.name))
  }
  else{
    groupIndex=messiDB$group.name==group
  }
  #### timeslices
  stopifnot(timeslice %in% c('all timeslices', timebins))
  if (timeslice == "all timeslices"){
    timesliceIndex=rep(TRUE,length(messiDB$Age))
  }
  else {
    timesliceIndex=messiDB$Age==timeslice
  }
  #### basisn
  stopifnot(basin %in% c("whole basin",regions))
  if (basin =="whole basin"){
    basinIndex=rep(TRUE,length(messiDB$region.new))
  }
  else {
    basinIndex=messiDB$region.new==basin
  }

  #### taxonomic level 
  stopifnot(taxLevel %in% c("species","genus","family"))
  if (taxLevel=="species"){
    taxIndex=!is.na(messiDB$Species.name) & !is.na(messiDB$Genus.name) & !is.na(messiDB$Family)
    occ=paste(messiDB$Genus.name, messiDB$Species.name, sep=' ')
    occ=occ[taxIndex ==TRUE & basinIndex==TRUE & timesliceIndex==TRUE & groupIndex==TRUE]
  }
  if (taxLevel=="genus"){
    taxIndex=!is.na(messiDB$Genus.name)  & !is.na(messiDB$Family)
    occ=paste(messiDB$Genus.name)
    occ=occ[taxIndex ==TRUE & basinIndex==TRUE & timesliceIndex==TRUE & groupIndex==TRUE]
  }
  
  if (taxLevel=="family"){
    taxIndex=!is.na(messiDB$Family)
    occ=paste(messiDB$Family)
    occ=occ[taxIndex ==TRUE & basinIndex==TRUE & timesliceIndex==TRUE & groupIndex==TRUE]
  }
  return(occ)
}

MSCSample(taxLevel = "species",
          group="all groups",
          timeslice = "all timeslices",
          basin = "whole basin")

rarefyTaxRichness=function(mySample, subsampleTo=1,noOfRep=1000){
  stopifnot(length(mySample)>=subsampleTo)
  stopifnot(!is.na(mySample))
  taxRichness=sapply(1:noOfRep, function(x) length(unique(sample(mySample,size=subsampleTo,replace=FALSE))))
  return(taxRichness)
}

rarefyEcoIndexes=function(mySample1,mySample2,subsampleTo=1,noOfRep=1000){
  stopifnot(length(mySample1)>=subsampleTo & length(mySample2)>=subsampleTo)
  stopifnot(!is.na(c(mySample1,mySample2)))
  
  out=list(soerensen=numeric(),simpson=numeric(),nestedness=numeric())
  for (i in 1:noOfRep){
    selectedocc1=sample(mySample1,size=subsampleTo,replace=FALSE)
    selectedocc2=sample(mySample2,size=subsampleTo,replace=FALSE)
    a=length(intersect(selectedocc1,selectedocc2))
    b=length(setdiff(selectedocc1,selectedocc2))
    c=length(setdiff(selectedocc2,selectedocc1))
    out$soerensen[i]=(b+c)/(2*a+b+c) # Baslega 2010
    out$simpson[i]=min(c(b,c))/(a+min(c(b,c)))
    out$nestedness[i]=(b+c)/(2*a+b+c)-(min(c(b,c))/(a+min(c(b,c))))
  }
  return(out)
}

}

#### eco indexes geographic comparison through time ####
## get sample isze
samplesize=c()
for (timebin in timebins){
  for (basin in regions){
    samplesize=c(samplesize,length(MSCSample(taxLevel='species',group='all groups',basin=basin,timeslice=timebin)))
  }
}
subsampleTo=floor(0.8*min(samplesize))
n=5000
for (timeslice in timebins){
  Eas=MSCSample(taxLevel='species',
                group="all groups",
                basin="Eastern Mediterranean",
                timeslice=timeslice)
  PoP=MSCSample(taxLevel='species',
                group="all groups",
                basin="Po Plain-Northern Adriatic",
                timeslice=timeslice)
  Wes=MSCSample(taxLevel='species',
                group="all groups",
                basin="Western Mediterranean",
                timeslice=timeslice)
  EasVsPoP=rarefyEcoIndexes(Eas,PoP,subsampleTo,noOfRep = n)
  PoPVsWes=rarefyEcoIndexes(PoP,Wes,subsampleTo,noOfRep = n)
  WesVsEas=rarefyEcoIndexes(Wes,Eas,subsampleTo,noOfRep = n)
  pdf(paste("Soerensen " ,timeslice ,", all groups, Species level, comparison between regions.pdf"))
  boxplot(list("East vs. Po Plain"=EasVsPoP$soerensen,
               "Po Plain vs. West"=PoPVsWes$soerensen,
               "West vs. East"=WesVsEas$soerensen),
          main=paste("Soerensen ", timeslice ,", all groups, Species level, comparison between regions"),
          ylab=paste("Subsampled to", subsampleTo, "occurrences"))
  dev.off()
  
  pdf(paste("Simpson " ,timeslice ,", all groups, Species level, comparison between regions.pdf"))
  boxplot(list("East vs. Po Plain"=EasVsPoP$simpson,
               "Po Plain vs. West"=PoPVsWes$simpson,
               "West vs. East"=WesVsEas$simpson),
          main=paste("Simpson ", timeslice ,", all groups, Species level, comparison between regions"),
          ylab=paste("Subsampled to", subsampleTo, "occurrences"))
  dev.off()
  
  pdf(paste("Nestedness " ,timeslice ,", all groups, Species level, comparison between regions.pdf"))
  boxplot(list("East vs. Po Plain"=EasVsPoP$nestedness,
               "Po Plain vs. West"=PoPVsWes$nestedness,
               "West vs. East"=WesVsEas$nestedness),
          main=paste("Nestedness ", timeslice ,", all groups, Species level, comparison between regions"),
          ylab=paste("Subsampled to", subsampleTo, "occurrences"))
  dev.off()
  print(timeslice)
}


#### species richenss
regions=c("Eastern Mediterranean","Po Plain-Northern Adriatic","Western Mediterranean")
for (timeslice in timebins){
  Eas=MSCSample(taxLevel='species',
                group="all groups",
                basin="Eastern Mediterranean",
                timeslice=timeslice)
  PoP=MSCSample(taxLevel='species',
                group="all groups",
                basin="Po Plain-Northern Adriatic",
                timeslice=timeslice)
  Wes=MSCSample(taxLevel='species',
                group="all groups",
                basin="Western Mediterranean",
                timeslice=timeslice)
  subsampleTo=900 #floor(max(1,0.8*min(c(length(Eas),length(PoP),length(Wes)))))
  print(subsampleTo)
  speciesRichnessList=list("Eastern Mediterranean"=rarefyTaxRichness(Eas,
                                                                     subsampleTo = subsampleTo,
                                                                     noOfRep = 50000),
                           "Po Plain-Northern Adriatic"=rarefyTaxRichness(PoP,
                                                                          subsampleTo = subsampleTo,
                                                                          noOfRep = 50000),
                           "Western Mediterranean"=rarefyTaxRichness(Wes,
                                                                     subsampleTo = subsampleTo,
                                                                     noOfRep = 50000))
  ymax=800 #max(sapply(speciesRichnessList,function(x) max(x)))
  pdf(paste("Species_richness " ,timeslice ,", all groups, Species level.pdf"))
  boxplot(speciesRichnessList,
          main=paste("Species Richness ", timeslice ,", all groups, Species level"),
          ylab=paste("Subsampled to", subsampleTo, "occurrences"),
          ylim=c(0,1.1*ymax))
  dev.off()
  
}


rarefyTaxRichness(a,subsampleTo = 10,noOfRep = 1000)

group.names=sort(unique(messiDB$group.name),na.last = TRUE)

#### fill in missing genus names with family names
has.genus.name=!is.na(messiDB$Genus.name)
messiDB$Genus.name[!has.genus.name]=paste(messiDB$Family[!has.genus.name]," indet.", sep='')

all(!is.na(messiDB$Genus.name))
n=10000

messiDB$Species.name[!is.na(messiDB$Species.name)]= paste(messiDB$Genus.name[!is.na(messiDB$Species.name)], messiDB$Species.name[!is.na(messiDB$Species.name)], sep=' ')


dataTo=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$Age==timebins[1]]
dataPre=messiDB$Species.name[!is.na(messiDB$Species.name)  & messiDB$Age==timebins[2]]
dataZa=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$Age==timebins[3]]

minN=min(length(dataTo),length(dataPre),length(dataZa))
sampleSize=max(1,floor(0.8*minN))
print(sampleSize)

TOME=myrarefaction(dataTo,dataPre,sampleSize,noofrep = n)
MEZA=myrarefaction(dataPre,dataZa,sampleSize,noofrep = n)

ymax=max(c(TOME$sr1,TOME$sr2,MEZA$sr2))
pdf(paste("Species_richness (Whole Basin), all groups, Species level.pdf"))
boxplot(list("Tortonian"=TOME$sr1,
             "pre-evaporitic Messinian"=TOME$sr2,
             "Zanclean"=MEZA$sr2 ),
        main=paste("Species Richness (Whole Basin) , all groups, Species level"),
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1.1*ymax))
dev.off()

pdf(paste("Soerensen (Whole Basin)  all groups Species level.pdf"))
boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$soerensen,
             "pre-evo Messinian\nvs.\nZanclean"=MEZA$soerensen),
        main="Soerensen (Whole Basin)   all groups  Species level",
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1))
dev.off()

pdf(paste("Simpson (Whole Basin)   all groups   Species level.pdf"))
boxplot(list("Tortonian\nvs.\npre-evo Messinian"=TOME$simpson,
             "pre-evo Messinian\nvs.\n Zanclean"=MEZA$simpson),
        main="Simpson (Whole Basin)   all groups   Species level",
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1))
dev.off()

pdf(paste("Nestedness (Whole Basin)  all groups  species level .pdf"))
boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$nestedness,
             "pre-evo Messinian\nvs.\n Zanclean"=MEZA$nestedness),
        main="Nestedness (Whole Basin)   all groups  Species level",
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1))
dev.off()

#### genus level ####
dataTo=messiDB$Genus.name[!is.na(messiDB$Genus.name) & messiDB$Age==timebins[1]]
dataPre=messiDB$Genus.name[!is.na(messiDB$Genus.name)  & messiDB$Age==timebins[2]]
dataZa=messiDB$Genus.name[!is.na(messiDB$Genus.name) & messiDB$Age==timebins[3]]

minN=min(length(dataTo),length(dataPre),length(dataZa))
sampleSize=max(1,floor(0.8*minN))
print(sampleSize)

TOME=myrarefaction(dataTo,dataPre,sampleSize,noofrep = n)
MEZA=myrarefaction(dataPre,dataZa,sampleSize,noofrep = n)

ymax=max(c(TOME$sr1,TOME$sr2,MEZA$sr2))
pdf(paste("Genus_richness (Whole Basin), all groups, Genus level.pdf"))
boxplot(list("Tortonian"=TOME$sr1,
             "pre-evaporitic Messinian"=TOME$sr2,
             "Zanclean"=MEZA$sr2 ),
        main=paste("Genus Richness (Whole Basin) , all groups, Genus level"),
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1.1*ymax))
dev.off()


#####

for ( group.name in group.names){
  # whole basin #
    dataTo=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$group.name == group.name & messiDB$Age==timebins[1]]
    dataPre=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$group.name == group.name & messiDB$Age==timebins[2]]
    dataZa=messiDB$Species.name[!is.na(messiDB$Species.name) & messiDB$group.name == group.name & messiDB$Age==timebins[3]]
    minN=min(length(dataTo),length(dataPre),length(dataZa))
    sampleSize=max(1,floor(0.8*minN))
    print(sampleSize)
    
    TOME=myrarefaction(dataTo,dataPre,sampleSize,noofrep = n)
    MEZA=myrarefaction(dataPre,dataZa,sampleSize,noofrep = n)
    
    ymax=max(c(TOME$sr1,TOME$sr2,MEZA$sr2))
    pdf(paste("Species_richness (Whole Basin)", group.name," Species level.pdf"))
    boxplot(list("Tortonian"=TOME$sr1,
                 "pre-evaporitic Messinian"=TOME$sr2,
                 "Zanclean"=MEZA$sr2 ),
            main=paste("Species Richness (Whole Basin) ", group.name," Species level"),
            ylab=paste("Subsampled to", sampleSize, "occurrences"),
            ylim=c(0,1.1*ymax))
    dev.off()
    
    pdf(paste("Soerensen (Whole Basin)  ", group.name," Species level.pdf"))
    boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$soerensen,
                 "pre-evo Messinian\nvs.\nZanclean"=MEZA$soerensen),
            main="Soerensen (Whole Basin)  ", group.name," Species level",
            ylab=paste("Subsampled to", sampleSize, "occurrences"),
            ylim=c(0,1))
    dev.off()
    
    pdf(paste("Simpson (Whole Basin)  ", group.name,"  Species level.pdf"))
    boxplot(list("Tortonian\nvs.\npre-evo Messinian"=TOME$simpson,
                 "pre-evo Messinian\nvs.\n Zanclean"=MEZA$simpson),
            main="Simpson (Whole Basin)  ", group.name,"  Species level",
            ylab=paste("Subsampled to", sampleSize, "occurrences"),
            ylim=c(0,1))
    dev.off()
    
    pdf(paste("Nestedness (Whole Basin) ", group.name,"species level .pdf"))
    boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$nestedness,
                 "pre-evo Messinian\nvs.\n Zanclean"=MEZA$nestedness),
            main="Nestedness (Whole Basin)  ", group.name," Species level",
            ylab=paste("Subsampled to", sampleSize, "occurrences"),
            ylim=c(0,1))
    dev.off()
    
    print(group.name)
}

tax.level=c("Family", "Genus.name","Species.name")
rd=list()
names=c("fish")
for(i in 1:length(names)){
  occ.name=paste(names[i],"_occ",sep='')
  loc.name=paste(names[i],"_loc",sep='')
  rd[[names[i]]]$occ=read.csv(paste(occ.name,".csv",sep=''),stringsAsFactors=FALSE)
  rd[[names[i]]]$loc=read.csv(paste(loc.name,".csv",sep=''),stringsAsFactors=FALSE)
  rd$fish$occ$Geographic.distribution=sapply(rd$fish$occ$Locality, function(x) rd$fish$loc$Region[x==rd$fish$loc$Locality])
}

#regional, sorted by times slices
rare_res1=list()
for (timebin in timebins){
  rare_res1[[timebin]]=list()
  for (tax in tax.level){
    dataE=rd$fish$occ[tax][rd$fish$occ$Age==timebin & rd$fish$occ$Geographic.distribution=="Eastern Mediterranean",]
    dataW=rd$fish$occ[tax][rd$fish$occ$Age==timebin & rd$fish$occ$Geographic.distribution=="Western Mediterranean",]
    dataP=rd$fish$occ[tax][rd$fish$occ$Age==timebin & rd$fish$occ$Geographic.distribution=="Po Plain-Northern Adriatic",]
    sampleSize=max(1,floor(0.8*min(length(dataE),length(dataW),length(dataP))))

    EW=myrarefaction(dataE,dataW,sampleSize,noofrep=1000)
    WP=myrarefaction(dataW,dataP,sampleSize,noofrep=1000)
    PE=myrarefaction(dataP,dataE,sampleSize,noofrep=1000)
    #rare_res1[timebin][tax]=list()
    rare_res1[[timebin]][[tax]]$EW=EW
    rare_res1[[timebin]][[tax]]$WP=WP
    rare_res1[[timebin]][[tax]]$PE=PE
    rare_res1[[timebin]][[tax]]$nosamples=sampleSize
    jpeg(paste("00fish species richness",tax, timebin, "comparison regions.jpeg"))
    boxplot(list(Eastern=rare_res1[[timebin]][[tax]]$EW$sr1,
                 Western=rare_res1[[timebin]][[tax]]$WP$sr1,
                 Poplain=rare_res1[[timebin]][[tax]]$PE$sr1),
            main=paste("species richness",tax, timebin))
    dev.off()
    
    jpeg(paste("00fish soerensen",tax, timebin, "comparison regions.jpeg"))
    boxplot(list(EvsW=rare_res1[[timebin]][[tax]]$EW$soerensen,
                 WvsP=rare_res1[[timebin]][[tax]]$WP$soerensen,
                 PvsE=rare_res1[[timebin]][[tax]]$PE$soerensen),
            main=paste("Soerensen",tax, timebin))
    dev.off()
    
    jpeg(paste("00fish simpson",tax, timebin, "comparison regions.jpeg"))
    boxplot(list(EvsW=rare_res1[[timebin]][[tax]]$EW$simpson,
                 WvsP=rare_res1[[timebin]][[tax]]$WP$simpson,
                 PvsE=rare_res1[[timebin]][[tax]]$PE$simpson),
            main=paste("Simpson",tax, timebin))
    dev.off()
    
    jpeg(paste("00fish nestedness",tax, timebin, "comparison regions.jpeg"))
    boxplot(list(EvsW=rare_res1[[timebin]][[tax]]$EW$nestedness,
                 WvsP=rare_res1[[timebin]][[tax]]$WP$nestedness,
                 PvsE=rare_res1[[timebin]][[tax]]$PE$nestedness),
            main=paste("Nestedness",tax, timebin))
    dev.off()
  }
}






# throug time, sorted by regions
rare_res2=list()
for (region in regions){
  rare_res2[[region]]=list()
  for (tax in tax.level){
    dataZ=rd$fish$occ[tax][rd$fish$occ$Age=="Zanclean" & rd$fish$occ$Geographic.distribution==region,]
    dataM=rd$fish$occ[tax][rd$fish$occ$Age=="pre-evaporitic Messinian" & rd$fish$occ$Geographic.distribution==region,]
    dataT=rd$fish$occ[tax][rd$fish$occ$Age=="Tortonian" & rd$fish$occ$Geographic.distribution==region,]

    sampleSize=max(1,floor(0.8*min(length(dataZ),length(dataM),length(dataT))))
    
    ZM=myrarefaction(dataZ,dataM,sampleSize,noofrep=1000)
    MT=myrarefaction(dataM,dataT,sampleSize,noofrep=1000)
    rare_res2[[region]][[tax]]$ZM=ZM
    rare_res2[[region]][[tax]]$MT=MT
    
    jpeg(paste("00fish species richness",tax, region, "comparison time.jpeg"))
    boxplot(list(Tortonian=rare_res2[[region]][[tax]]$MT$sr2,
                 peMessinian=rare_res2[[region]][[tax]]$MT$sr1,
                 Zanclean=rare_res2[[region]][[tax]]$ZM$sr1 ),
            main=paste("species richness",tax, region))
    dev.off()
    
    jpeg(paste("00fish soerensen",tax, region, "comparison time.jpeg"))
    boxplot(list(TortonianvsppeMessinian=rare_res2[[region]][[tax]]$MT$soerensen,
                 peMessinianvsZanclean=rare_res2[[region]][[tax]]$ZM$soerensen ),
            main=paste("soerensen",tax, region))
    dev.off()
    
    jpeg(paste("00fish simpson",tax, region, "comparison time.jpeg"))
    boxplot(list(TortonianvsppeMessinian=rare_res2[[region]][[tax]]$MT$simpson,
                 peMessinianvsZanclean=rare_res2[[region]][[tax]]$ZM$simpson ),
            main=paste("simpson",tax, region))
    dev.off()
    
    jpeg(paste("00fish nestedness",tax, region, "comparison time.jpeg"))
    boxplot(list(TortonianvsppeMessinian=rare_res2[[region]][[tax]]$MT$nestedness,
                 peMessinianvsZanclean=rare_res2[[region]][[tax]]$ZM$nestedness ),
            main=paste("nestedness",tax, region))
    dev.off()
    
  }
}


















