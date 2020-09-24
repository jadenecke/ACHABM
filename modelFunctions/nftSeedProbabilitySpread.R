##########################################################
# 
##########################################################

nftSeedProbabilitySpread <- function(neurons, lag, maxSeedProb, seedProbMaxRise, nftName = "nft", nftSeedProbabilityName = "nftSeedProbability", nftSeedProbabilityStackName = "nftSeedProbabilityStack", connectionsName = "conMat"){
  require(parallel)
  spreadList <- nftSeedProbabilitySpread_calculateTransfer(neurons, lag, nftName)
  index <- sapply(spreadList, function(l) {sum(l[[1]]) > 0})
  if(sum(index) > 1){
    neurons <- nftSeedProbabilitySpread_transferFromNeuronToNeuron(neurons, spreadList, index, nftSeedProbabilityStackName, connectionsName)  
  }
  for(i in seq_along(neurons)){
    neurons[[i]][[nftSeedProbabilityName]] <- min(neurons[[i]][[nftSeedProbabilityName]] + min(neurons[[i]][[nftSeedProbabilityStackName]], seedProbMaxRise), maxSeedProb)
    neurons[[i]][[nftSeedProbabilityStackName]] <- 0
  }
  return(neurons)
}


nftSeedProbabilitySpread_calculateTransfer <- function(neurons, lag, nftName){
  return(parallel::mclapply(neurons, function(l){
    if(l[["alive"]]){
    spreadVec <- rbeta(length(l[[nftName]]),
                       sqrt(l[[nftName]]),
                        lag) * (l[["activity"]] != 0)
    # if(any(is.na(spreadVec))){l; spreadVec}
    return(spreadVec)
    } else {
      return(c(0,0))
    }
  })) #, mc.preschedule = TRUE, mc.cores = 2
}

nftSeedProbabilitySpread_transferFromNeuronToNeuron <- function(neurons, spreadList, index, nftSeedProbabilityStackName, connectionsName){
  for(i in seq(length(spreadList[index]))){
    con <- neurons[index][[i]][[connectionsName]]
    for(j in seq_along(con)){
      neurons[[con[j]]][[nftSeedProbabilityStackName]] <- neurons[[con[j]]][[nftSeedProbabilityStackName]] + spreadList[index][[i]][j]
    }
  }
  return(neurons)
}



# curve(.7/(1 + exp(-4 * (.5 + x))), 0, 1)
