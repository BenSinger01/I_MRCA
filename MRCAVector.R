#class of vectors storing locations of MRCAs of pairs of leaves
setClass("MRCAVector", representation(LeafNo = "numeric", MRCALocations = "character"))
setMethod("initialize", "MRCAVector", function(.Object, ..., MRCALocations){callNextMethod(.Object, ..., MRCALocations = MRCALocations, LeafNo =  ((1+sqrt(1+8*length(MRCALocations)))/2))})

#class of vectors storing nodes of MRCAs of pairs of leaves, and probabilities over possible locations of the nodes
setClass("ProbabilisticMRCAVector", representation(LeafNo = "numeric", NodeLocationProbabilities = "matrix", MRCANodes = "character"))
setMethod("initialize", "ProbabilisticMRCAVector", function(.Object, ..., MRCANodes, NodeLocationProbabilities){callNextMethod(.Object, ..., NodeLocationProbabilities = NodeLocationProbabilities, MRCANodes = MRCANodes, LeafNo =  ((1+sqrt(1+8*length(MRCANodes)))/2))})

#class of geographic networks, storing vertices and distances
setClass("Network", representation(Locations = "character", Distances = "matrix"))

#function to calculate incompatibility of two MRCAVectors
I_MRCA <- function(v1,v2,N) {
  #check that the vectors are of the same class
  if(class(v1) != class(v2)){
    stop("Vectors are of different classes - check if one is probabilistic and the other not")
  }
  #if the vectors are not probabilistic, just sum the distances on the network
  #of corresponding elements in the MRCA vectors
  else if(class(v1) == "MRCAVector"){
    I <- 0
    for (i in 1:length(v1@MRCALocations)){
      I <- I + N@Distances[v1@MRCALocations[i],v2@MRCALocations[i]]
    }
    return(I)
  }
  #if the vectors are probabilistic, sum distances weighted by probabilities,
  #as outlined in the notes/paper accompanying this code
  else if(class(v1) == "ProbabilisticMRCAVector"){
    I <- 0
    locPairs = combn(N@Locations,2)
    for (i in 1:length(v1@MRCANodes)){
      #sum over pairs of different vertices (don't need to sum over identical
      #vertices since each vertex as a distance of zero from itself)
      for (j in 1:(length(locPairs)/2)){
        I <- I + v1@NodeLocationProbabilities[v1@MRCANodes[i],locPairs[1,j]]*v2@NodeLocationProbabilities[v2@MRCANodes[i],locPairs[2,j]]*N@Distances[locPairs[1,j],locPairs[2,j]]
        I <- I + v1@NodeLocationProbabilities[v1@MRCANodes[i],locPairs[2,j]]*v2@NodeLocationProbabilities[v2@MRCANodes[i],locPairs[1,j]]*N@Distances[locPairs[1,j],locPairs[2,j]]
      }
    }
    return(I)
  }
  #return an error message if the inputs are the wrong type
  else {
    stop("Vectors are not of the class MRCAVector or ProbabilisticMRCAVector")
  }
}

#Tips is here to make the Get*Vector functions a bit prettier
Tips <- function(beast){
  return(unlist(beast@phylo["tip.label"]))
}

#GetMRCAVector: function to extract a deterministic MRCA location vector from
#a data frame generated using the read.beast function. The comparison argument
#can be assigned to another treedata object in order to make sure that only
#shared leaves are used in creating the vector
GetMRCAVector <- function(beast, comparison = 0){
  #if there's a comparison, get rid of the tips not shared between the trees
  if (class(comparison) == "treedata"){
    leaves <- Tips(beast)[Tips(beast) %in% Tips(comparison)]
  }
  else{
    leaves <- Tips(beast)
  }
  #create a vector of all possible pairs of leaf names
  pairs <- combn(leaves,2)
  n_choose <- length(pairs)/2
  #initialise the location vector
  V <- character(n_choose)
  #for each pair of leaves, assign the corresponding entry in V to the location
  #of those leaves' MRCA in the BEAST-generated tree
  for (i in 1:n_choose){
    V[i] <- toString(beast@data[MRCA(beast,pairs[,i]),"Location"])
  }
  #return a deterministic MRCA vector
  return(new("MRCAVector", MRCALocations=V))
}

#GetProbabilisticMRCAVecotor: function to extract a probabilistic MRCA location
#vector from a data frame generated using the read.beast function
GetProbabilisticMRCAVector <- function(beast, comparison = 0){
  #if there's a comparison, get rid of the tips not shared between the trees
  if (class(comparison) == "treedata"){
    leaves <- Tips(beast)[Tips(beast) %in% Tips(comparison)]
    locations <- unique(c(unlist(beast@data["Location.set"]),unlist(comparison@data["Location.set"])))
  }
  else{
    leaves <- Tips(beast)
    locations <- unique(unlist(beast@data["Location.set"]))
  }
  #find some useful quantities
  n_leaves <- length(leaves)
  n_nodes <- n_leaves-1
  #create a vector of all possible pairs of leaf names
  pairs <- combn(leaves,2)
  n_choose <- length(pairs)/2
  #initialise the node vector
  V <- character(n_choose)  #for each pair of leaves, assign the corresponding entry in V to the index
  #of the MRCA node
  for (i in 1:n_choose){
    V[i] <- MRCA(beast,pairs[,i])
  }
  #find the set of nodes
  nodes = unlist(beast@data["node"])[unlist(beast@data["node"]) %in% V]
  #initialise the location probability array
  X <- array(0,dim = c(n_nodes,length(locations)), dimnames = list(nodes,locations))
  #for each node, enter the corresponding location entries into X
  for (n in nodes){
    locs <- unlist(beast@data[n,"Location.set"])
    for (j in 1:length(locs)){
      X[n,locs[j]] <- unlist(beast@data[n,"Location.set.prob"])[j]
    }
  }
  return(new("ProbabilisticMRCAVector", NodeLocationProbabilities = X, MRCANodes = V))
}

#GetNetwork: returns a Network object from a data frame made with read.beast,
#with all locations distance 1 away from each other
GetNework <- function(beast, comparison = 0){
  if (class(comparison) == "treedata"){
    locations <- unique(c(unlist(beast@data["Location.set"]),unlist(comparison@data["Location.set"])))
  }
  else{
    ocations <- unique(unlist(beast@data["Location.set"]))
  }
  distances <- 1-diag(length(locations))
  rownames(distances) <- locations
  colnames(distances) <- locations
  return(new("Network", Locations = locations, Distances = distances))
}
