##############################################
# Generate a dummy network with numberOfCluster circular areas connected with each other as adjacency list.
# 28.08.2020
# By: Jannis Denecke

genNClusterNetwork <- function(numberOfCluster, clusterMinSize, clusterMaxSize, dendritesPerNeuronMin, dendritesPerNeuronMax, synapsesPerDendriteMin, synapsesPerDendriteMax){
  # Paramter Checks
  if(synapsesPerDendriteMax * dendritesPerNeuronMax > clusterMinSize) stop("Possibly to many connections: Either reduce synapsesPerDendriteMin//dendritesPerNeuronMin and/or increase ClusterMinSize")
  require(parallel)
  require(rlist)
  # Calculate Coordinates for Neurons
  (size <- ceiling(runif(numberOfCluster, clusterMinSize, clusterMaxSize))) # number of neurons in each area
  # connectionMatrix <- matrix(floor(runif(numberOfCluster^2,
  #                                        min(size) * synapsesPerDendriteMin * dendritesPerNeuronMin,
  #                                        max(size) * synapsesPerDendriteMax * dendritesPerNeuronMax)),
  #                            nrow = numberOfCluster)
  connectionMatrix <- matrix(0, nrow = numberOfCluster, ncol= numberOfCluster)
  for(i in seq_along(size)){
    for(j in seq_along(size)){
      if(i != j){
        connectionMatrix[i,j] <- runif(1,
                                       min(c(size[i], size[j])) * synapsesPerDendriteMin * dendritesPerNeuronMin,
                                       min(c(size[i], size[j])) * synapsesPerDendriteMax * dendritesPerNeuronMax
                                       )
      }
    }
  }
  
  
  (connectionMatrix <- floor((connectionMatrix + t(connectionMatrix)) / 2)) #symmetrification of the matrix
  # connectionMatrix[,1] <- connectionMatrix[,1] /4
  
  
  coords <- expand.grid(1:ceiling((clusterMaxSize*2.2)^(1/3)),
                        1:ceiling((clusterMaxSize*2.2)^(1/3)),
                        1:ceiling((clusterMaxSize*2.2)^(1/3))
  )
  coords <- coords - ceiling(ceiling((clusterMaxSize*2.2)^(1/3))/2)
  coords <- coords[order(apply(coords, 1, function(x){sqrt(x %*% x)})), ] #order by euclidean distance
  coords <- as.matrix(coords)
  
  coordsArea <- expand.grid(1:(ceiling(numberOfCluster^(1/3))), 1:ceiling((numberOfCluster^(1/3))), 1:ceiling((numberOfCluster^(1/3))))
  coordsArea <- coordsArea - ceiling(ceiling(numberOfCluster^(1/3))/2)
  coordsArea <- coordsArea[order(apply(coordsArea, 1, function(x){sqrt(x %*% x)})), ]
  coordsArea <- as.matrix(coordsArea)
  
  # Visualize the possible coordinates of one area / all areas:
  # library("rgl")
  # lim <- 60000
  # plot3d(coords[1:lim, 1],
  #        coords[1:lim, 2],
  #        coords[1:lim, 3],
  #        aspect = FALSE)
  # plot3d(unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],1] + coordsArea[i,1] * ceiling((clusterMaxSize*2.2)^(1/3))})),
  #        unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],2] + coordsArea[i,2] * ceiling((clusterMaxSize*2.2)^(1/3))})),
  #        unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],3] + coordsArea[i,3] * ceiling((clusterMaxSize*2.2)^(1/3))})),
  #        aspect = FALSE)
  
  #initiate Model:
  neurons <- initiateModel(size = size,
                           connectionMatrix = connectionMatrix,
                           dendritesPerNeuronMin = dendritesPerNeuronMin,
                           dendritesPerNeuronMax = dendritesPerNeuronMax,
                           synapsesPerDendriteMin = synapsesPerDendriteMin,
                           synapsesPerDendriteMax = synapsesPerDendriteMax)
  
  #add single coordinates:
  neurons <- mapply(c, neurons, x = unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],1] + coordsArea[i,1] * ceiling((clusterMaxSize*2.2)^(1/3))})), SIMPLIFY = FALSE)
  neurons <- mapply(c, neurons, y = unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],2] + coordsArea[i,2] * ceiling((clusterMaxSize*2.2)^(1/3))})), SIMPLIFY = FALSE)
  neurons <- mapply(c, neurons, z = unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],3] + coordsArea[i,3] * ceiling((clusterMaxSize*2.2)^(1/3))})), SIMPLIFY = FALSE)
  #add coordinate vector:
  neurons <- mapply(c,
                    neurons,
                    coords = apply(cbind(unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],1] + coordsArea[i,1] * ceiling((clusterMaxSize*2.2)^(1/3))})),
                                         unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],2] + coordsArea[i,2] * ceiling((clusterMaxSize*2.2)^(1/3))})),
                                         unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],3] + coordsArea[i,3] * ceiling((clusterMaxSize*2.2)^(1/3))}))),
                                   1,
                                   list),
                    SIMPLIFY = FALSE)
  
  
  
  
}



######################################### deprecated ######################################

##### Network / Igraph approach:
# g <- graph.empty(sum(size), directed = FALSE)
# g <- set_vertex_attr(g, "area", value = paste0("a", rep(1:numberOfCluster, size)))
# g <- set_vertex_attr(g, "x", value = unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],1] + coordsArea[i,1] * ceiling((clusterMaxSize*2.2)^(1/3))})))
# g <- set_vertex_attr(g, "y", value = unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],2] + coordsArea[i,2] * ceiling((clusterMaxSize*2.2)^(1/3))})))
# g <- set_vertex_attr(g, "z", value = unlist(lapply(1:numberOfCluster, function(i){coords[1:size[i],3] + coordsArea[i,3] * ceiling((clusterMaxSize*2.2)^(1/3))})))

# for(i in 1:(numberOfCluster-1)){
#   for(j in (i+1):numberOfCluster){
#     print(paste0(i, "-", j))
#     print(connectionMatrix[i,j])
#     connectionsPerVertex <- floor(connectionMatrix[i,j] / size[i])
#     originVertices <- as.vector(V(g)[get.vertex.attribute(g, "area") == paste0("a", i)])
#     targetVertices <- as.vector(V(g)[get.vertex.attribute(g, "area") == paste0("a", j)])
#     
#     newEdges <- sapply(originVertices,
#                        function(v){
#                          as.vector(t(cbind(v, sample(targetVertices, connectionsPerVertex))))
#                        }
#     )
#     g <- add_edges(g, newEdges)
#   }
# }

# write.graph(g, "4Cluster.graphml", format = "graphml")
# plot(g, 
#      layout = t(apply(cbind(get.vertex.attribute(g, "x"),
#                             get.vertex.attribute(g, "y"),
#                             get.vertex.attribute(g, "z")
#      ), 1, function(v) rotateCoordinatesByX(rotateCoordinatesByZ(rotateCoordinatesByY(v, 30), 30), 60))),
#      vertex.size = 1,
#      vertex.label = NA,
#      shape = "circle",
#      edge.label = NA,
#      edge.color = rgb(.1,.1,.1, .1)
# )

# as.vector(t(matrix(1:6, ncol= 2)))

