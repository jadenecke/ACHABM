initiateModel <- function(size, connectionMatrix, dendritesPerNeuronMin, dendritesPerNeuronMax, synapsesPerDendriteMin, synapsesPerDendriteMax){
  require(parallel)
  require(rlist)
  #Create empty List of Lists with one list for each neuron (id = 1:numberOfNeurons)
  neurons <- lapply(1:sum(size), function(x){list(id = x)})
  names(neurons) <- 1:sum(size)
  neurons <- mapply(c, neurons, area = rep(1:length(size), size), SIMPLIFY = FALSE) #add area id to each neuron
  if(!is.null(names(size))){
    neurons <- mapply(c, neurons, areaName = rep(names(size), size), SIMPLIFY = FALSE)
  } else {
    neurons <- mapply(c, neurons, areaName = rep(paste0("area_", seq(length(size))), size), SIMPLIFY = FALSE)
  }
  
  
  conList <- lapply(1:sum(size), function(x){list()})
  #add connection Matrix:
  for(i in 1:length(size)){
    print(paste0("Connecting areas: ", i))
    targetAreas <- setdiff(1:length(size), i)
    print(paste0("Connecting to target areas: ", paste0(targetAreas,  collapse = " / ")))
    print(paste0("Area size: ", size[i], " / Target area sizes: ", paste0(size[targetAreas], collapse = " / ")))
    
    connectionsPerNeuronIntended <- connectionMatrix[i, targetAreas] / size[i]
    originVertices <- if(i == 1){
      (1:size[1])
    } else {
      (sum(size[1:(i-1)])+1) : sum(size[1:i])
    }
    
    
    
    for(j in targetAreas){
      print(j)
      originIndex <- floor(runif(connectionMatrix[i, j], min = min(originVertices), max = max(originVertices) +1))
      if(j == 1){
        targetVertices <- (1:size[1])
      } else {
        targetVertices <- (sum(size[1:(j-1)])+1) : sum(size[1:j])
      }
      
      for(k in seq_along(originIndex)){
        conList[[originIndex[k]]][[paste0("conArea", j)]] <- list.append(conList[[originIndex[k]]][[paste0("conArea", j)]], sample(targetVertices, 1))
      }
    }
  }
  
  cMat <- outer(dendritesPerNeuronMin:dendritesPerNeuronMax, synapsesPerDendriteMin:synapsesPerDendriteMax, "*")
  dendriteNumberList <- pbmcapply::pbmclapply(conList, function(l) {
    len <- sapply(l, length)
    cBest <- matrix(NA, length(len), dim(cMat)[2])
    for(i in seq_along(len)){
      cBest[i,] <- apply(len[i] - cMat, 2, function(x) withCallingHandlers(suppressWarnings(min(x[x >= 0]))))
    }
    return((synapsesPerDendriteMin:synapsesPerDendriteMax)[which.min(colSums(cBest, na.rm = TRUE))])
  })
  
  pb = txtProgressBar(min = 0, max = length(neurons), initial = 0, style = 3) 
  for(i in seq_along(neurons)){
    setTxtProgressBar(pb,i)
    neurons[[i]][["conMat"]] <- matrix(unlist(lapply(conList[[i]], function(v){
      l <- floor(length(v) / dendriteNumberList[[i]])
      return(matrix(v[seq(l * dendriteNumberList[[i]])], ncol = dendriteNumberList[[i]], nrow = l))
    })), ncol = dendriteNumberList[[i]], byrow = T)
    neurons[[i]][["dendriteAreas"]] <- rep(as.numeric(gsub(".*?([0-9]+).*", "\\1", names(conList[[i]]))),
                                           unlist(lapply(conList[[i]], function(v){floor(length(v) / dendriteNumberList[[i]])}))
                                           )
  }
  
  
  #Count number of connections for each area:
  conMatReal <- matrix(0, dim(connectionMatrix)[1], dim(connectionMatrix)[2])
  pb = txtProgressBar(min = 0, max = length(neurons), initial = 0, style = 3) 
  for(i in seq_along(neurons)){
    setTxtProgressBar(pb,i)
    v <- table(neurons[[i]][["dendriteAreas"]])  * dim(neurons[[i]][["conMat"]])[2]
    o <- neurons[[i]][["area"]]
    d <- as.numeric(names(v))
    for(j in seq_along(d)){
      conMatReal[o, d[j]] <- conMatReal[o, d[j]] + v[j]
    }
    
  }
  print("Connection Matrix reduced to fit model:")
  print(conMatReal)
  print("Percent reduction:")
  print(round((1 - (conMatReal/connectionMatrix)) * 100, 1))
  
  
  
  #####################################
  
  # for(i in 1:length(size)){
  #   print(paste0("Connecting areas: ", i))
  #   targetAreas <- setdiff(1:length(size), i)
  #   print(paste0("Connecting to target areas: ", paste0(targetAreas,  collapse = " / ")))
  #   print(paste0("Area size: ", size[i], " / Target area sizes: ", paste0(size[targetAreas], collapse = " / ")))
  #   
  #   connectionsPerNeuronIntended <- connectionMatrix[i, targetAreas] / size[i]
  #   originVertices <- if(i == 1){
  #     (1:size[1])
  #   } else {
  #     (sum(size[1:(i-1)])+1) : sum(size[1:i])
  #   }
  #   
  #   # dpnMinMax <- connectionsPerNeuronIntended / dendritesPerNeuronMin
  #   # spdMin <- max(ceiling(connectionsPerNeuronIntended / dendritesPerNeuronMin))
  #   # synapsesPerDendrite <- spdMin + which.min(lapply(spdMin:synapsesPerDendriteMax, function(s){sum(abs(dpnMinMax %% s))})) - 1
  #   # dendritesPerNeuron <- floor(connectionsPerNeuronIntended / synapsesPerDendrite)
  #   # connectionsPerNeuron <- synapsesPerDendrite * dendritesPerNeuron
  #   
  #   #print(paste0("New outbound conections for this area (deviation): ", paste((connectionsPerNeuron * size[i]) - connectionMatrix[i,targetAreas], collapse = " / ")))
  #   print(paste0("New outbound conections for this area (intended): ", paste0(connectionMatrix[i,targetAreas], collapse = " / ")))
  #   # print(paste0("New outbound conections for this area (deviation absolute): ", paste0(connectionsPerNeuron * size[i] - connectionMatrix[i,targetAreas], collapse = " / ")))
  #   # print(paste0("New outbound conections for this area (deviation percent): ", paste0("-", round(1 - (connectionsPerNeuron * size[i] / connectionMatrix[i,targetAreas]), 3) * 100, "%", collapse = " / ")))
  #   # print(paste0("New outbound conections for per Neuron (intended): ", paste0(connectionsPerNeuronIntended, collapse = " / ")))
  #   # print(paste0("New outbound conections for per Neuron (deviation): ", paste0(synapsesPerDendrite * dendritesPerNeuron  - connectionsPerNeuronIntended, collapse = " / ")))
  #   # print(paste0("Actual new outbound connections (dendrites per Neuron): ", paste0(dendritesPerNeuron, collapse = " / ")))
  #   # print(paste0("Actual new outbound connections (snypases per dendrites): ", paste0(synapsesPerDendrite, collapse = " / ")))
  #   
  #   
  #   
  #   
  #   
  #   conListPerArea <- lapply(seq_along(targetAreas), function(ti){
  #     if(targetAreas[ti] == 1){
  #       targetVertices <- (1:size[1])
  #     } else {
  #       targetVertices <- (sum(size[1:(targetAreas[ti]-1)])+1) : sum(size[1:targetAreas[ti]])
  #     }
  #     print(ti)
  #     cl <- lapply(originVertices, function(o) sample(targetVertices, connectionsPerNeuron[ti]))
  #     return(matrix(unlist(cl), nrow = length(originVertices), byrow = T))
  #   })
  #   
  #   conListPerNeuron <- lapply(split(do.call(cbind, conListPerArea), seq(size[i])),
  #                              function(v){matrix(v, ncol = synapsesPerDendrite, nrow = sum(dendritesPerNeuron), byrow = TRUE)})
  #   
  #   neurons[originVertices] <- mapply(list.append,
  #                                     neurons[originVertices],
  #                                     conMat = conListPerNeuron, SIMPLIFY = FALSE)
  # }
  
  
  #Add Amyloid and NFT load Matrices
  templateList <- lapply(neurons, function(l){
    d <- dim(l$conMat)
    matrix(integer(1), nrow = d[1], ncol = d[2])
  })
  
  # t1 <- lapply(neurons, function(l){
  #   d <- dim(l$conMat)
  #   array(integer(1), c(d, 3))
  # })
  
  #somewhat safe additions to the framework
  neurons <- mapply(list.append, neurons, aMonomer = templateList, SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, aDimer = templateList, SIMPLIFY = FALSE)
  #neurons <- mapply(list.append, neurons, aAggregate = templateList, SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, aPlaque = templateList, SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, nft = templateList, SIMPLIFY = FALSE)
  
  #current additions (might be removed):
  neurons <- mapply(list.append, neurons, activity = lapply(templateList, `+`, 1), SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, aMonomerStack = as.list(rep(0,length(neurons))), SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, aDimerStack = as.list(rep(0,length(neurons))), SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, aAggregateCount = templateList, SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, aAggregateSum = templateList, SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, nftSeedProbability = as.list(rep(nftGeneration_nftSeedProbability,length(neurons))), SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, nftSeedProbabilityStack = as.list(rep(0,length(neurons))), SIMPLIFY = FALSE)
  neurons <- mapply(list.append, neurons, alive = as.list(rep(TRUE,length(neurons))), SIMPLIFY = FALSE)
  
  return(neurons)
}
