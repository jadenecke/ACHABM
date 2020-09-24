##########################################################
# 
##########################################################

synpaseActivityUpdate <- function(neurons,
                                  declineStability,
                                  declineStart,
                                  activityCutOff,
                                  nftCutOffMean,
                                  aliveName = "alive",
                                  nftName = "nft",
                                  activityName = "activity",
                                  aAggregateCountName = "aAggregateCount",
                                  aDimerName =  "aDimer"
){
  return(mclapply(neurons, function(l){
    if(l[["alive"]]){
      activity <- rbeta(length(l[[activityName]]),
                        declineStability / sqrt((l[[aDimerName]] + l[[aAggregateCountName]] + l[[nftName]])),
                        (l[[aDimerName]] + l[[aAggregateCountName]] + l[[nftName]])^2 / declineStart
      )
      
      activity[activity <= activityCutOff] <- 0
      activity[l[[nftName]] >= nftCutOffMean] <- 0
      l[[activityName]] <- matrix(activity, dim(l[[activityName]])[1], dim(l[[activityName]])[2])
      
      if(mean(l[[activityName]]) <= activityCutOff | mean(l[[nftName]]) > nftCutOffMean){
        l[[aliveName]] <- FALSE
        #l[[activityName]] <- matrix(0, dim(l[[activityName]])[1], dim(l[[activityName]])[2])
      }
    }
    return(l)
  }))
  
}

#### Way to add activity reduction by neighbouring plaques: Not sure whether i want to add this right now, because i think this could be already a result from the existing rules
# m <- matrix(0, 7, 10)
# m[4, 5] <- 100
# 
#
# #### BAND MATRIX:
# plaquePullIntraDendrite <- 4
# if(ceiling(plaquePullIntraDendrite / 2) >= (dim(m)[2] - 1)){
#   mt <- matrix(1, nrow = dim(m)[2], ncol = dim(m)[2])
# } else {
#   mt <- diag(1, nrow = dim(m)[2], ncol = dim(m)[2])
#   for(i in seq(ceiling(plaquePullIntraDendrite / 2))){
#     diag(mt[0:(-i),]) <- 1
#     diag(mt[,0:(-i)]) <- 1
#   }
# }
# if(ceiling(plaquePullIntraDendrite / 2) >= (dim(m)[1] - 1)){
#   mt2 <- matrix(1, nrow = dim(m)[1], ncol = dim(m)[1])
# } else {
#   mt2 <- diag(1, nrow = dim(m)[1], ncol = dim(m)[1])
#   for(i in seq(ceiling(plaquePullIntraDendrite / 2))){
#     diag(mt2[0:(-i),]) <- 1/i
#     diag(mt2[,0:(-i)]) <- 1/i
#   }
# }
# t(t(m %*% mt) %*% mt2)

