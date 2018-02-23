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
        print(locPairs[2,j])
        I <- I + v1@NodeLocationProbabilities[v1@MRCANodes[i],locPairs[1,j]]*v2@NodeLocationProbabilities[v2@MRCANodes[i],locPairs[2,j]]*N@Distances[locPairs[1,j],locPairs[2,j]]
      }
    }
    return(I)
  }
  #return an error message if the inputs are the wrong type
  else {
    stop("Vectors are not of the class MRCAVector or ProbabilisticMRCAVector")
  }
}

#GetMRCAVecotor: function to extract a deterministic MRCA location vector from
#a data frame generated using the read.beast function
GetMRCAVector <- function(beast){
  #assign useful quantities
  n_nodes <- dim(beast@data)[1]
  n_leaves <- (n_nodes+1)/2
  #create a vector of all possible pairs of leaf indices
  pairs <- combn(1:n_leaves,2)
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
GetProbabilisticMRCAVector <- function(beast){
  #assign useful quantities
  n_nodes <- dim(beast@data)[1]
  n_leaves <- (n_nodes+1)/2
  locations <- unique(unlist(beast@data["Location.set"]))
  #create a vector of all possible pairs of leaf indices
  pairs <- combn(1:n_leaves,2)
  n_choose <- length(pairs)/2
  #initialise the node vector
  V <- character(n_choose)
  #initialise the location probability array
  X <- array(0,dim = c(n_nodes,length(locations)), dimnames = list(1:n_nodes,locations))
  #for each pair of leaves, assign the corresponding entry in V to the index
  #of the MRCA node
  for (i in 1:n_choose){
    V[i] <- MRCA(beast,pairs[,i])
  }
  #for each node, enter the corresponding location entries into X
  for (i in 1:n_nodes){
    locs <- unlist(beast@data[i,"Location.set"])
    for (j in 1:length(locs)){
      X[i,locs[j]] <- unlist(beast@data[i,"Location.set.prob"])[j]
    }
  }
  return(new("ProbabilisticMRCAVector", NodeLocationProbabilities = X, MRCANodes = V))
}

#GetNetwork: returns a Network object from a data frame made with read.beast,
#with all locations distance 1 away from each other
GetNework <- function(beast){
  locations <- unique(unlist(beast@data["Location.set"]))
  distances <- 1-diag(length(locations))
  rownames(distances) <- locations
  colnames(distances) <- locations
  return(new("Network", Locations = locations, Distances = distances))
}