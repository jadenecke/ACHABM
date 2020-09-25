##########################################################
# 
##########################################################

logModel <- function(neurons, logPath, step){
  #copy Model parameter
  if(step == 1){
    file.copy("modelParameter.R", file.path(logPath, "modelParameter.R"))
  }
  
  #start logging
  log <- list(step = step,
              neurons = mclapply(neurons, function(l){
                return(list(
                  id = l$id,
                  area = l$area,
                  areaName = l$areaName,
                  coords = ifelse(is.null(l$coords), no = l$coords, yes = NA),
                  numberOfsynapses = length(l$conMat),
                  connectionsTo = l$dendriteAreas,
                  aMonomerSum = sum(l$aMonomer),
                  aDimerSum = sum(l$aDimer),
                  aAggregateCountSum = sum(l$aAggregateCount),
                  aAggregateSumSum = sum(l$AggregateSum),
                  aPlaqueNumber = sum(l$aPlaque != 0),
                  aPlaqueSum = sum(l$aPlaque),
                  nftSeedProbability = l$nftSeedProbability,
                  nftSum = sum(l$nft),
                  alive = l$alive,
                  activity = mean(l$activity),
                  synapseAlivePercent = mean(l$activity != 0)
                ))
              })
              
  )
  
  saveRDS(log, file = file.path(logPath, paste0("modelLog_", step, ".rds")))
  
  return(NULL)
}
