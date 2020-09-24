##########################################################
# 
##########################################################

nftGeneration <- function(neurons, nftAcceleration, probabilitySaclingFactor, maxNftGrowth, synapseActivityUpdate_nftCutOffMean, nftFlatClearance, nftName = "nft", aAggregateCountName = "aAggregateCount", nftSeedProbability = "nftSeedProbability", activityName = "activity", aDimerName = "aDimer"){
  return(mclapply(neurons, function(l){
    if(l[["alive"]]){
      l[[nftName]] <- pmin(l[[nftName]] + rbinom(length(l[[nftName]]),
                                                 pmin(ceiling(sqrt(l[[aAggregateCountName]] + l[[nftName]] + l[[aDimerName]]) / nftAcceleration),
                                                      (synapseActivityUpdate_nftCutOffMean * maxNftGrowth)),
                                                 l[[nftSeedProbability]]
      ),
      2^31,
      na.rm = TRUE
      )
      # l[[nftName]] <- pmin(l[[nftName]] + rbinom(length(l[[nftName]]),
      #                                            ceiling(sqrt(l[[aAggregateCountName]] + l[[nftName]] + l[[aDimerName]]) / nftAcceleration),
      #                                            l[[nftSeedProbability]] * probabilitySaclingFactor
      # ),
      # 2^31,
      # na.rm = TRUE
      # )
    #NFT CLEANING:  
      l[[nftName]] <- pmax(l[[nftName]] - rbinom(length(l[[nftName]]),
                                                 nftFlatClearance,
                                                 .5
      ),
      0,
      na.rm = TRUE
      )
    }
    return(l)
  }))

}


# plot(unlist(lapply(seq(1, 10000, 1), function(x){rbinom(1,
#                                                           pmin(ceiling(sqrt(x) / nftAcceleration), (20000 * .002)),
#                                                           l[[nftSeedProbability]] * probabilitySaclingFactor
# )})))
