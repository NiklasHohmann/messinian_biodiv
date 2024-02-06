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
