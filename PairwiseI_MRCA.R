Beast <- function(name){
  return (read.beast(file = gsub("y",name,"/homes/singer/Documents/FMDV Data/MCCtree_n74_y.txt")))
}

PairwiseI_MRCA <- function(names, option = "none"){
  size = length(names)
  dists = array(rep(0,2*(size**2)),dim = c(size,size,2),dimnames = list(names,names,c("I","n.mutual.leaves")))
  if (option == "normalise"){ x = TRUE } else { x = FALSE }
  
  for (i in 1:(length(names)-1)){
    for (j in (1+i):length(names)){
      beast1 = Beast(names[i])
      beast2 = Beast(names[j])
      v1 = GetProbMRCAVector(beast1,comparison=beast2)
      v2 = GetProbMRCAVector(beast2,comparison=beast1)
      N = GetNetwork(beast1,comparison=beast2)
      leaf.logic = Tips(beast1) %in% Tips(beast2)
      dists[i,j,"n.mutual.leaves"] = sum(leaf.logic)
      dists[i,j,"I"] = I_MRCA(v1,v2,N,normalisation=x)
    }
  }
  return(dists)
}