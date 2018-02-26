# This code is to calculate the MRCA-matching distance between two phylogeographies.
# It follows the Google Style Guide for R.

# Classes

# class of vectors storing locations of MRCAs of pairs of leaves
setClass("MRCA.vector", representation(leaf.no = "numeric", MRCA.location = "character"))
setMethod("initialize", "MRCA.vector", function(.Object, ..., MRCA.location)
  {callNextMethod(.Object, ..., MRCA.location = MRCA.location, 
                  leaf.no =  ((1+sqrt(1+8*length(MRCA.location)))/2))})

# class of vectors storing nodes of MRCAs of pairs of leaves, and probabilities over possible locations of the nodes
setClass("prob.MRCA.vector", 
         representation(leaf.no = "numeric", node.location.prob = "matrix", 
                        MRCA.node = "character"))
setMethod("initialize", "prob.MRCA.vector", 
          function(.Object, ..., MRCA.node, node.location.prob)
            {callNextMethod(.Object, ..., 
                            node.location.prob = node.location.prob, MRCA.node = MRCA.node,
                            leaf.no =  ((1+sqrt(1+8*length(MRCA.node)))/2))})

# class of geographic networks, storing vertices and distances
setClass("network", representation(locations = "character", distances = "matrix"))

# Functions

# Tips is here to make later functions more readable
Tips <- function(beast){
  return(unlist(beast@phylo["tip.label"]))
}

# MRCAVector: function to extract a deterministic MRCA location vector from
# a data frame generated using the read.beast function. The comparison argument
# can be assigned to another treedata object in order to make sure that only
# shared leaves are used in creating the vector
MRCAVector <- function(beast, comparison = 0){
  #i f there's a comparison, get rid of the tips not shared between the trees
  if (class(comparison) == "treedata"){
    leaves <- Tips(beast)[Tips(beast) %in% Tips(comparison)]
  }
  else{
    leaves <- Tips(beast)
  }
  # create a vector of all possible pairs of leaf names
  pairs <- combn(leaves,2)
  nChoose <- length(pairs)/2
  # initialise the location vector
  V <- character(nChoose)
  # for each pair of leaves, assign the corresponding entry in V to the location
  # of those leaves' MRCA in the BEAST-generated tree
  for (i in 1:nChoose){
    V[i] <- toString(beast@data[MRCA(beast,pairs[,i]),"Location"])
  }
  # return a deterministic MRCA vector
  return(new("MRCA.vector", MRCA.location=V))
}

# ProbMRCAVector: function to extract a probabilistic MRCA location
# vector from a data frame generated using the read.beast function
ProbMRCAVector <- function(beast, comparison = 0){
  # if there's a comparison, get rid of the tips not shared between the trees
  if (class(comparison) == "treedata"){
    leaves <- Tips(beast)[Tips(beast) %in% Tips(comparison)]
    locs <- unique(c(unlist(beast@data["Location.set"]),
                          unlist(comparison@data["Location.set"])))
  }
  else{
    leaves <- Tips(beast)
    locs <- unique(unlist(beast@data["Location.set"]))
  }
  # find some useful quantities
  nLeaves <- length(leaves)
  nNodes <- nLeaves-1
  # create a vector of all possible pairs of leaf names
  pairs <- combn(leaves,2)
  nChoose <- length(pairs)/2
  # initialise the node vector
  V <- character(nChoose)  
  # for each pair of leaves, assign the corresponding entry in V to the index
  # of the MRCA node
  for (i in 1:nChoose){
    V[i] <- MRCA(beast,pairs[,i])
  }
  # find the set of nodes
  nodes = unlist(beast@data["node"])[unlist(beast@data["node"]) %in% V]
  # initialise the location probability array
  X <- array(0,dim = c(nNodes,length(locs)), 
             dimnames = list(nodes,locs))
  # for each node, enter the corresponding location entries into X
  for (n in nodes){
    l <- unlist(beast@data[n,"Location.set"])
    for (j in 1:length(l)){
      X[n,l[j]] <- unlist(beast@data[n,"Location.set.prob"])[j]
    }
  }
  return(new("prob.MRCA.vector", node.location.prob = X,
             MRCA.node = V))
}

# Network: returns a network object from a data frame made with read.beast,
# with all locations distance 1 away from each other
Network <- function(beast, comparison = 0){
  if (class(comparison) == "treedata"){
    locs <- unique(c(unlist(beast@data["Location.set"]),
                          unlist(comparison@data["Location.set"])))
  }
  else{
    locs <- unique(unlist(beast@data["Location.set"]))
  }
  dists <- 1-diag(length(locs))
  rownames(dists) <- locs
  colnames(dists) <- locs
  return(new("network", locations = locs, distances = dists))
}

# function to calculate incompatibility of two MRCA.vectors
I_MRCA <- function(v1,v2,N) {
  # check that the vectors are of the same class
  if(class(v1) != class(v2)){
    stop("Vectors are of different classes - check if one is probabilistic and the other not")
  } else if(class(v1) == "MRCA.vector") {
    # if the vectors are not probabilistic, just sum the distances on the network
    # of corresponding elements in the MRCA vectors
    I <- 0
    for (i in 1:length(v1@MRCA.location)){
      I <- I + N@distances[v1@MRCA.location[i],v2@MRCA.location[i]]
    }
    return(I)
  } else if(class(v1) == "prob.MRCA.vector") {
    # if the vectors are probabilistic, sum distances weighted by probabilities,
    # as outlined in the notes/paper accompanying this code
    I <- 0
    locPairs = combn(N@locations,2)
    for (i in 1:length(v1@MRCA.node)){
      # sum over pairs of different vertices (don't need to sum over identical
      # vertices since each vertex as a distance of zero from itself). Sum
      # performed for each order of each pair.
      for (j in 1:(length(locPairs)/2)){
        I <- I + v1@node.location.prob[v1@MRCA.node[i],locPairs[1,j]]*
          v2@node.location.prob[v2@MRCA.node[i],locPairs[2,j]]*
          N@distances[locPairs[1,j],locPairs[2,j]]
        I <- I + v1@node.location.prob[v1@MRCA.node[i],locPairs[2,j]]*
          v2@node.location.prob[v2@MRCA.node[i],locPairs[1,j]]*
          N@distances[locPairs[1,j],locPairs[2,j]]
      }
    }
    return(I)
  } else {
    # return an error message if the inputs are the wrong type
    stop("Vectors are not of the class MRCA.vector or ProbabilisticMRCA.vector")
  }
}