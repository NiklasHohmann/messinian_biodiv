# Maps of localities - Agiadi et al. "A revised marine fossil record of the Mediterranean before and after the Messinian Salinity Crisis"

# 21.02.2024
# Author: K. Agiadi

library(tidyverse)
library(ggmap)
library(ggthemes)
library(ggplot2)


# To get the Mediterranean map from Stadia Maps, it is necessary to have an API key, which you can get for free by registering at https://client.stadiamaps.com/
# Then, register the key by
# register_stadiamaps("PUT YOUR API KEY HERE", write=TRUE)
# To setup the bbox margins, go to OpenStreetMap, find the view you want and click export

bbox=c(left=-7.4,bottom=28.5,right=40,top=47.1)
Med_map<-get_stadiamap(bbox,maptype="stamen_toner_background",zoom=5)  
ggmap(Med_map)

# Read locality data
coord<-read.csv("coord.csv")
class(coord)

# Whole dataset map
qmplot(lon,lat,data=coord,maptype="stamen_toner_lite",color=I("red"))
# OR
ggmap(Med_map)+geom_point(data=coord,aes(color="darkred",shape=Age))
ggsave("figs/maps/dataset_map.pdf",width=9,height=5) 

# Color by stage
Tortonian<-filter(coord,Age=="Tortonian")
Messinian<-filter(coord,Age=="pre-evaporitic Messinian")
Zanclean<-filter(coord,Age=="Zanclean")
ggmap(Med_map)+geom_point(data=Tortonian,color="darkblue")+geom_point(data=Messinian,color="darkred")+geom_point(data=Zanclean,color="gold2")
ggsave("figs/maps/locs_by_stage_map.pdf",width=9,height=5)

# Color by region
wMed<-filter(coord,region=="Western Mediterranean")
eMed<-filter(coord,region=="Eastern Mediterranean")
PoA<-filter(coord,region=="Po Plain-Northern Adriatic")
ggmap(Med_map)+geom_point(data=wMed,color="green4")+geom_point(data=eMed,color="pink4")+geom_point(data=PoA,color="orange3")
ggsave("figs/maps/locs_by_region_map.pdf",width=9,height=5)

# Color by stage and region
ggmap(Med_map)+
  geom_point(data=wMed,color="turquoise4",aes(shape=Age))+
  geom_point(data=eMed,color="firebrick2",aes(shape=Age))+
  geom_point(data=PoA,color="gold1",aes(shape=Age))+
  labs(x="longitude",y="latitude")+
  theme(legend.position="bottom")

ggsave("figs/maps/stage_region_map2.pdf",width=9,height=5)

# plot for each group separately
tax_groups<-unique(coord$group.name)
pforams<-filter(coord,group.name=="planktic_foraminifera")
nano<-filter(coord,group.name=="nanoplankton")
ostracods<-filter(coord,group.name=="ostracods")
bforams<-filter(coord,group.name=="benthic_foraminifera")
gastr<-filter(coord,group.name=="gastropods")
echinoids<-filter(coord,group.name=="echinoids")
fish<-filter(coord,group.name=="fish")
mmammals<-filter(coord, group.name=="marine_mammals")
biv<-filter(coord,group.name=="bivalves")
corals<-filter(coord,group.name=="corals")
dino<-filter(coord,group.name=="dinocysts")
sharks<-filter(coord,group.name=="sharks")
bryo<-filter(coord,group.name=="bryozoans")
othermolluscs<-filter(coord, group.name=="scaphopod_chitons_cephalopods")

ggmap(Med_map)+geom_point(data=pforams,aes(color=Age))
ggsave("figs/maps/pforams.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=nano,aes(color=Age))
ggsave("figs/maps/nano.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=ostracods,aes(color=Age))
ggsave("figs/maps/ostracods.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=bforams,aes(color=Age))
ggsave("figs/maps/bforams.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=gastr,aes(color=Age))
ggsave("figs/maps/gastro.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=echinoids,aes(color=Age))
ggsave("figs/maps/echinoids.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=fish,aes(color=Age))
ggsave("figs/maps/fish.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=mmammals,aes(color=Age))
ggsave("figs/maps/marine_mammals.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=biv,aes(color=Age))
ggsave("figs/maps/bivalves.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=corals,aes(color=Age))
ggsave("figs/maps/corals.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=dino,aes(color=Age))
ggsave("figs/maps/dino.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=sharks,aes(color=Age))
ggsave("figs/maps/sharks.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=bryo,aes(color=Age))
ggsave("figs/maps/bryo.pdf",width=9,height=5)

ggmap(Med_map)+geom_point(data=othermolluscs,aes(color=Age))
ggsave("figs/maps/othermolluscs.pdf",width=9,height=5)

# plot groups split into plankton, benthos, necton and large marine vertebrates
# plankton: planktic_foraminifera, nanoplankton, dinocysts
# benthos: benthic_foraminifera, ostracods, bryozoans, molluscs, echinoids, corals
# necton: fish 
# large marine vertebrates: sharks, marine mammals
plankton<-rbind(pforams,nano,dino)
benthos<-rbind(bforams,ostracods,bryo,othermolluscs,gastr,biv,echinoids,corals)
necton<-rbind(fish)
large_mar_vert<-rbind(sharks,mmammals)

ggmap(Med_map)+geom_point(data=plankton,color="gold1",aes(shape=Age))+geom_point(data=benthos,color="red",aes(shape=Age))+geom_point(data=necton,color="blue",aes(shape=Age))+geom_point(data=large_mar_vert,color="green",aes(shape=Age))
ggsave("figs/maps/fgroups.pdf",width=9,height=5)
