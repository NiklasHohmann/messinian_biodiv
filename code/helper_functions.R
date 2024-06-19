get_from_db=function(group, basin, timeslice, taxLevel = "species", tax_groups = "joint corals"){
  #'
  #' @title extract occurrence from Messinian database
  #' 
  #' @description
  #' returns a vector of taxon names from the Messinian DB on the specified taxonomic level, region, and time slice
  #' 
  #' 
  #' @param taxLevel "species","genus" or "family". Taxonomic level to extract
  #' @param group character, element of _group.names_ or "all groups". Taxonomic groups to retreive
  #' @param basin character, element of _regions_ or "whole basin". From which area occurrences are selected
  #' @param timeslice character, element of _timebins_ or "all timeslices". Time interval of interest
  #' @param tax_groups type of taxonomic groups used. either "joint corals" or "split corals". determines whether corals are split into a and z corals or not
  #' 
  #' @returns character vector of taxon names on the specified taxonomic level
  
  #### groups 
  stopifnot(tax_groups %in% c("joint corals", "split corals"))
  if (tax_groups == "joint corals"){
    stopifnot(group %in% group.names.ext)
  }
  if (tax_groups == "split corals"){
    stopifnot(group %in% group.names.ext_sc)
  }
  
  if (group=='all groups'){
    groupIndex=rep(TRUE,length(messinian_db$group.name))
  }
  else{
    if (tax_groups == "joint corals"){
      groupIndex=messinian_db$group.name==group
    }
    if (tax_groups == "split corals"){
      groupIndex = messinian_db$group.name_sc == group
    }
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
    occ=messinian_db$Genus.name
    occ=occ[taxIndex ==TRUE & basinIndex==TRUE & timesliceIndex==TRUE & groupIndex==TRUE]
  }
  
  if (taxLevel=="family"){
    taxIndex=!is.na(messinian_db$Family)
    occ=messinian_db$Family
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

rarefyTaxGradient=function(mySample1, mySample2, subsampleTo,noOfRep){
  #'
  #' @title rarefy differences in taxonomic richness
  #' 
  #' @description
  #' use rarefaction do determine the difference in taxonomic richness between two samples
  #' 
  #' @param mySample1 vector of taxonomic names in first sample
  #' @param mySample2 vector of taxon names in second sample
  #' @param subsampleTo integer, nomber of occurrences to target for subsampling. Must be larger than length(mySample)
  #' @param noOfRep integer, number of subsampling repetitions
  #' 
  #' @returns integer vector of length noOfRep, containing tax richnesses at each repetition
  stopifnot(length(mySample1)>=subsampleTo)
  stopifnot(length(mySample2)>=subsampleTo)
  stopifnot(!is.na(mySample1))
  stopifnot(!is.na(mySample2))
  gradient=sapply(seq_len(noOfRep), function(x) length(unique(sample(mySample1,size=subsampleTo,replace=FALSE))) - length(unique(sample(mySample2,size=subsampleTo,replace=FALSE))))
  return(gradient)
}
