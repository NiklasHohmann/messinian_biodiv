get_from_db=function(group, basin, timeslice, taxLevel = "species"){
  #'
  #' @title extract occurrence from Messinian database
  #' 
  #' @description
  #' returns a vector of taxon names from the Messinian DB on the specified taxonomic level, sub basin, and time slice
  #' 
  #' 
  #' @param taxLevel "species","genus" or "family". Taxonomic level to extract
  #' @param group character, element of _group.names_ or "all groups". Taxonomic groups to retreive
  #' @param basin character, element of _regions_ or "whole basin". From which area occurrences are selected
  #' @param timeslice character, element of _timebins_ or "all timeslices". Time interval of interest
  #' 
  #' @returns character vector of taxon names on the specified taxonomic level
  
  #### groups 
  stopifnot(group %in% c('all groups',group.names))
  if (group=='all groups'){
    groupIndex=rep(TRUE,length(messinian_db$group.name))
  }
  else{
    groupIndex=messinian_db$group.name==group
  }
  #### timeslices
  stopifnot(timeslice %in% c('all timeslices', timebins))
  if (timeslice == "all timeslices"){
    timesliceIndex=rep(TRUE,length(messinian_db$Age))
  }
  else {
    timesliceIndex=messinian_db$Age==timeslice
  }
  #### basisn
  stopifnot(basin %in% c("whole basin",regions))
  if (basin =="whole basin"){
    basinIndex=rep(TRUE,length(messinian_db$region.new))
  }
  else {
    basinIndex=messinian_db$region.new==basin
  }
  
  #### taxonomic level 
  stopifnot(taxLevel %in% c("species","genus","family"))
  if (taxLevel=="species"){
    taxIndex=!is.na(messinian_db$Species.name) & !is.na(messinian_db$Genus.name) & !is.na(messinian_db$Family)
    occ=paste(messinian_db$Genus.name, messinian_db$Species.name, sep=' ')
    occ=occ[taxIndex ==TRUE & basinIndex==TRUE & timesliceIndex==TRUE & groupIndex==TRUE]
  }
  if (taxLevel=="genus"){
    taxIndex=!is.na(messinian_db$Genus.name)  & !is.na(messinian_db$Family)
    occ=paste(messinian_db$Genus.name)
    occ=occ[taxIndex ==TRUE & basinIndex==TRUE & timesliceIndex==TRUE & groupIndex==TRUE]
  }
  
  if (taxLevel=="family"){
    taxIndex=!is.na(messinian_db$Family)
    occ=paste(messinian_db$Family)
    occ=occ[taxIndex ==TRUE & basinIndex==TRUE & timesliceIndex==TRUE & groupIndex==TRUE]
  }
  return(occ)
}


rarefyTaxRichness=function(mySample, subsampleTo,noOfRep){
  #'
  #' @title rarefy taxonomic richness
  #' 
  #' @description
    #' rarefies a vector of taxon names to a specific sample size
    #' 
  #' @param mySample vector of taxonomic names
  #' @param subsampleTo integer, nomber of occurrences to target for subsampling. Must be larger than length(mySample)
  #' @param noOfRep integer, number of subsampling repetitions
  #' 
  #' @returns integer vector of length noOfRep, containing tax richnesses at each repetition
  stopifnot(length(mySample)>=subsampleTo)
  stopifnot(!is.na(mySample))
  taxRichness=sapply(seq_len(noOfRep), function(x) length(unique(sample(mySample,size=subsampleTo,replace=FALSE))))
  return(taxRichness)
}

rarefyEcoIndexes=function(mySample1, mySample2, subsampleTo, noOfRep){
  #' 
  #' @title pairwise rarefaction for ecological parameters
  #' 
  #' @param mySample1 first vector of taxon names
  #' @param mySample2 second vector of taxon names
  #' @param subsampleto target sample size for subsampling.
  #' @param noOfRep integer, number of subsampling repetitions
  #' 
  #' @description
  #' performs pairwise subsampling from two vectors of taxon names, and returns 
  #' key ecological indices (soerensen & simpson index, nestedness)
  #'
  #' @returns a list with three names elements: "soerensen", "simpson", "nestedness", each a vector of length _noOfRep_, containing the ecological indices of the i-th subsampling run 
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
























biodiv_ind = function(v1, v2, n, no_of_rep){
  #' @title get biodiversity indices
  #' 
  #' @param v1 first vector of taxon names
  #' @param v2 second vector of taxon names
  #' @param n samplesize for rarefaction
  #' @param no_of_rep number of repetitions for rarefaction
  #' 
  #' @returns a list with elements names "tr1" and "tr2" (taxonomic richness of v1 and v2), "soerensen", "simpson", and "nestedness"
  #' 
  out=list("tr1"=numeric(length = n),
           "tr2"=numeric(length = n),
           "soerensen"=numeric(length = n),
           "simpson"=numeric(length = n),
           "nestedness"=numeric(length = n))
  
  for (i in seq_len(no_of_rep)){
    a = joint_rarefaction(v1,v2,n)
    x1 = a[[1]]
    x2 = a[[2]]
    out[["tr1"]][i] = tax_richness(x1)
    out[["tr2"]][i] = tax_richness(x2)
    out[["soerensen"]][i] = soerensen(x1,x2)
    out[["simpson"]][i] = simpson(x1,x2)
    out[["nestedness"]][i] = nestedness(x1,x2) 
  }
  return(out)
}

joint_rarefaction = function(v1, v2, n){
  #' @title rarefy from 2 vectors simultaneously
  #' 
  #' @description
  #' randomly selects n elements from 2 vectors in parallel 
  #' 
  #' 
  #' @param v1 first vector
  #' @param v2 second vector
  #' @param n no of elements drawn from vectors
  #' 
  #' @return a list with two elements, each is a vector of length n drawn from v1 or v2
  x1 = sample(x = v1, size = n, replace = FALSE)
  x2 = sample(x = v2, size = n, replace = FALSE)
  li = list(x1, x2)
  return(li)
}

nestedness = function(x1, x2){
  #' @title nestedness
  #'
  #' @param x1 vector of taxon names
  #' @param x2 vector of taxon names
  #' 
  #' @returns numeric, nestedness
  a=length(intersect(x1,x2))
  b=length(setdiff(x1,x2))
  c=length(setdiff(x2,x1))
  
  nestedness = (b+c)/(2*a+b+c)-(min(c(b,c))/(a+min(c(b,c))))
  return(nestedness)
}

simpson = function(x1, x2){
  #' @title simpson index
  #' 
  #' @param x1 vector of taxon names
  #' @param x2 vector of taxon names
  #' 
  #' @returns numeric, simpson index
  
  a=length(intersect(x1,x2))
  b=length(setdiff(x1,x2))
  c=length(setdiff(x2,x1))
  
  simpson_index=min(c(b,c))/(a+min(c(b,c)))
  return(simpson_index)
}

soerensen = function(x1,x2){
  #' @title soerensen index
  #' 
  #' @param x1 vector of taxon names
  #' @param x2 vector of taxon names
  #' 
  #' @returns numeric, soerensen index
  
  a=length(intersect(x1,x2))
  b=length(setdiff(x1,x2))
  c=length(setdiff(x2,x1))
  soerensen_index = (b+c)/(2*a+b+c)
  return(soerensen__index)
}

tax_richness = function(x){
  #' @title taxonomic richness
  #' 
  #' @param x vector of taxon names
  #' 
  #' @returns integer, taxonomic richness
  tax_richness = length(unique(x))
  return(tax_richness)
}
