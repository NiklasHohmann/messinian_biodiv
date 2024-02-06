#### Load data ####
messinian_db <- read.csv(file = "data/messinianDB.csv")

mdb = data.frame(id = messinian_db$ID,group = messinian_db$group.name)

#### Define constants ####
time_bins <- unique(messinian_db$Age)[c(3, 1, 2)] # sorted from old to young
regions <- unique(messinian_db$region.new)
group_names <- unique(messinian_db$group.name)


#### subsampling

## what species id status??
# what I need:
# sr all groups individually & combined - each interval
# biodiv measures for all groups individually & combined - comparing all pairs of intervals

# sr for all groups individually & combined for each subbasin
# biodiv measures for all groups individually & combined - comparing all pairs of intervals for each subbasin

for (group in group_names){
  Tor  = messinian_db$Species.name[messinian_db$group.name == group & messinian_db$Age == "Tortonian" ]
  Mes = messinian_db$Species.name[messinian_db$group.name == group & messinian_db$Age == "pre-evaporitic Messinian" ]
  Zan = messinian_db$Species.name[messinian_db$group.name == group & messinian_db$Age == "Zanclean" ]

}


#### fill in missing genus names with family names
has.genus.name=!is.na(messiDB$Genus.name)
messiDB$Genus.name[!has.genus.name]=paste(messiDB$Family[!has.genus.name]," indet.", sep='')

all(!is.na(messiDB$Genus.name))
n=10000

for ( group.name in group.names){
  # whole basin #
    dataTo=messiDB$Genus.name[messiDB$group.name == group.name & messiDB$Age==timebins[1]]
    dataPre=messiDB$Genus.name[messiDB$group.name == group.name & messiDB$Age==timebins[2]]
    dataZa=messiDB$Genus.name[messiDB$group.name == group.name & messiDB$Age==timebins[3]]
    minN=min(length(dataTo),length(dataPre),length(dataZa))
    sampleSize=floor(0.8*minN)
    print(sampleSize)
    
    TOME=myrarefaction(dataTo,dataPre,sampleSize,noofrep = n)
    MEZA=myrarefaction(dataPre,dataZa,sampleSize,noofrep = n)
    
    ymax=max(c(TOME$sr1,TOME$sr2,MEZA$sr2))
    pdf(paste("Species_richness (Whole Basin)", group.name," Genus level.pdf"))
    boxplot(list("Tortonian"=TOME$sr1,
                 "pre-evaporitic Messinian"=TOME$sr2,
                 "Zanclean"=MEZA$sr2 ),
            main=paste("Species Richness (Whole Basin) ", group.name," Genus level"),
            ylab=paste("Subsampled to", sampleSize, "occurrences"),
            ylim=c(0,1.1*ymax))
    dev.off()
    
    pdf(paste("Soerensen (Whole Basin)  ", group.name," Genus level.pdf"))
    boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$soerensen,
                 "pre-evo Messinian\nvs.\nZanclean"=MEZA$soerensen),
            main="Soerensen (Whole Basin)  ", group.name," Genus level",
            ylab=paste("Subsampled to", sampleSize, "occurrences"),
            ylim=c(0,1))
    dev.off()
    
    pdf(paste("Simpson (Whole Basin)  ", group.name,"  Genus level.pdf"))
    boxplot(list("Tortonian\nvs.\npre-evo Messinian"=TOME$simpson,
                 "pre-evo Messinian\nvs.\n Zanclean"=MEZA$simpson),
            main="Simpson (Whole Basin)  ", group.name,"  Genus level",
            ylab=paste("Subsampled to", sampleSize, "occurrences"),
            ylim=c(0,1))
    dev.off()
    
    pdf(paste("Nestedness (Whole Basin) ", group.name," .pdf"))
    boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$nestedness,
                 "pre-evo Messinian\nvs.\n Zanclean"=MEZA$nestedness),
            main="Nestedness (Whole Basin)  ", group.name," Genus level",
            ylab=paste("Subsampled to", sampleSize, "occurrences"),
            ylim=c(0,1))
    dev.off()
    
    print(group.name)
}


#### all groups together ####

# whole basin #
dataTo=messiDB$Genus.name[ messiDB$Age==timebins[1]]
dataPre=messiDB$Genus.name[ messiDB$Age==timebins[2]]
dataZa=messiDB$Genus.name[ messiDB$Age==timebins[3]]
minN=min(length(dataTo),length(dataPre),length(dataZa))
sampleSize=floor(0.8*minN)
print(sampleSize)

TOME=myrarefaction(dataTo,dataPre,sampleSize,noofrep = n)
MEZA=myrarefaction(dataPre,dataZa,sampleSize,noofrep = n)

ymax=max(c(TOME$sr1,TOME$sr2,MEZA$sr2))
pdf(paste("Species_richness (Whole Basin) all groups Genus level.pdf"))
boxplot(list("Tortonian"=TOME$sr1,
             "pre-evaporitic Messinian"=TOME$sr2,
             "Zanclean"=MEZA$sr2 ),
        main=paste("Species Richness (Whole Basin) all groups Genus level"),
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1.1*ymax))
dev.off()

pdf(paste("Soerensen (Whole Basin)  all groups Genus level.pdf"))
boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$soerensen,
             "pre-evo Messinian\nvs.\nZanclean"=MEZA$soerensen),
        main="Soerensen (Whole Basin)  all groups Genus level",
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1))
dev.off()

pdf(paste("Simpson (Whole Basin)  all groups  Genus level.pdf"))
boxplot(list("Tortonian\nvs.\npre-evo Messinian"=TOME$simpson,
             "pre-evo Messinian\nvs.\n Zanclean"=MEZA$simpson),
        main="Simpson (Whole Basin)  all groups  Genus level",
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1))
dev.off()

pdf(paste("Nestedness (Whole Basin) all groups .pdf"))
boxplot(list("Tortonian \nvs.\npre-evo Messinian"=TOME$nestedness,
             "pre-evo Messinian\nvs.\n Zanclean"=MEZA$nestedness),
        main="Nestedness (Whole Basin)  all groups Genus level",
        ylab=paste("Subsampled to", sampleSize, "occurrences"),
        ylim=c(0,1))
dev.off()



#### pl ####

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
















myrarefaction=function(dataset1,dataset2,sampleSize,noofrep=1000){
  noofrep=max(1,floor(noofrep))
  out=list(sr1=numeric(length = ),sr2=numeric(),soerensen=numeric(),simpson=numeric(),nestedness=numeric())
  for (i in 1:noofrep){
    selectedocc1=sample(dataset1,size=sampleSize,replace=FALSE)
    selectedocc2=sample(dataset2,size=sampleSize,replace=FALSE)
    out$sr1[i]=length(unique(selectedocc1))
    out$sr2[i]=length(unique(selectedocc2))
    a=length(intersect(selectedocc1,selectedocc2))
    b=length(setdiff(selectedocc1,selectedocc2))
    c=length(setdiff(selectedocc2,selectedocc1))
    out$soerensen[i]=(b+c)/(2*a+b+c)
    out$simpson[i]=min(c(b,c))/(a+min(c(b,c)))
    out$nestedness[i]=(b+c)/(2*a+b+c)-(min(c(b,c))/(a+min(c(b,c))))
  }
  return(out)
}

