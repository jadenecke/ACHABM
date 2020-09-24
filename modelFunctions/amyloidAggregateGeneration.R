##########################################################
# 
##########################################################

amyloidAggregateGeneration <- function(neurons,
                                       seedingProbabiltyMax,
                                       seedingProbabilityCurveSteepness,
                                       seedingProbabilityCurveInflectionPoint,
                                       aggregateGrowthProbabiltyMax,
                                       aggregateGrowthProbabilityCurveSteepness,
                                       aggregateGrowthDelay,
                                       maxDecline,
                                       stability,
                                       aggregateMaxSize = 20,
                                       matMonomerName = "aMonomer",
                                       matDimerName = "aDimer",
                                       matAggregateCountName = "aAggregateCount",
                                       matAggregateSumName = "aAggregateSum"){
  require("parallel")
  index <- sapply(neurons, function(l) {sum(l[[matAggregateCountName]]) > 0})
  neurons[index] <- amyloidAggregateGeneration_aggregateGrowthByDimer(neurons[index], aggregateGrowthProbabiltyMax, aggregateGrowthProbabilityCurveSteepness, aggregateGrowthDelay, aggregateMaxSize, matDimerName, matAggregateCountName, matAggregateSumName)
  neurons[index] <- amyloidAggregateGeneration_aggregateGrowthByMonomer(neurons[index], aggregateGrowthProbabiltyMax, aggregateGrowthProbabilityCurveSteepness, aggregateGrowthDelay, aggregateMaxSize, matMonomerName, matAggregateCountName, matAggregateSumName)
  neurons[index] <- amyloidAggregateGeneration_aggregateDissaggregationIntoMonomer(neurons[index], maxDecline, stability, matMonomerName, matAggregateCountName, matAggregateSumName, aggregateMaxSize)
  #neurons <- amyloidAggregateGeneration_aggregateDeseeding(neurons, matMonomerName, matAggregateCountName, matAggregateSumName)
  neurons <- amyloidAggregateGeneration_aggregateSeeding(neurons, seedingProbabiltyMax, seedingProbabilityCurveSteepness, seedingProbabilityCurveInflectionPoint, matMonomerName,  matDimerName, matAggregateCountName, matAggregateSumName)
  return(neurons)
  }

amyloidAggregateGeneration_aggregateSeeding <- function(neurons, seedingProbabiltyMax, seedingProbabilityCurveSteepness, seedingProbabilityCurveInflectionPoint, matMonomerName,  matDimerName, matAggregateCountName, matAggregateSumName){
  return(mclapply(neurons, function(l){
    if(l[["alive"]]){
      newAgg <- rbinom(n = length(l[[matMonomerName]]),
                       size = pmin(l[[matMonomerName]], l[[matDimerName]]),
                       prob = seedingProbabiltyMax/(1 + exp(-seedingProbabilityCurveSteepness*((l[[matMonomerName]] + 2 * l[[matDimerName]]) - seedingProbabilityCurveInflectionPoint)))
      )
      if(sum(newAgg) != 0){
        l[[matAggregateCountName]] <- l[[matAggregateCountName]] + newAgg
        l[[matAggregateSumName]] <- l[[matAggregateSumName]] + (newAgg * 3)
        l[[matMonomerName]] <- l[[matMonomerName]] - newAgg
        l[[matDimerName]] <- l[[matDimerName]] - newAgg
      }
    }
    return(l)
  }))
}

amyloidAggregateGeneration_aggregateGrowthByMonomer <- function(neurons, aggregateGrowthProbabiltyMax, aggregateGrowthProbabilityCurveSteepness, aggregateGrowthDelay, aggregateMaxSize, matMonomerName, matAggregateCountName, matAggregateSumName){
  return(mclapply(neurons, function(l){
    openAggregates <- l[[matAggregateCountName]] * aggregateMaxSize - l[[matAggregateSumName]]
    transfer <- rbinom(n = length(l[[matAggregateSumName]]),
                       size = pmin(openAggregates, l[[matMonomerName]]),
                       prob = aggregateGrowthProbabiltyMax/(1 + exp(-aggregateGrowthProbabilityCurveSteepness*((l[[matMonomerName]]) - (aggregateGrowthDelay - (floor(openAggregates / aggregateMaxSize)^2)))))
    )
    l[[matAggregateSumName]] <- l[[matAggregateSumName]] + transfer
    l[[matMonomerName]] <- l[[matMonomerName]] - transfer
    return(l)
  }))
}

amyloidAggregateGeneration_aggregateGrowthByDimer <- function(neurons, aggregateGrowthProbabiltyMax, aggregateGrowthProbabilityCurveSteepness, aggregateGrowthDelay, aggregateMaxSize, matDimerName, matAggregateCountName, matAggregateSumName){
  return(mclapply(neurons, function(l){
    openAggregates <- floor((l[[matAggregateCountName]] * aggregateMaxSize - l[[matAggregateSumName]]) / 2)
    transfer <- rbinom(n = length(l[[matAggregateSumName]]),
                       size = pmin(openAggregates, l[[matDimerName]]),
                       prob = aggregateGrowthProbabiltyMax/(1 + exp(-aggregateGrowthProbabilityCurveSteepness*((l[[matDimerName]]) - (aggregateGrowthDelay - (floor(openAggregates / aggregateMaxSize)^2)))))
    )
    l[[matAggregateSumName]] <- l[[matAggregateSumName]] + 2 * transfer
    l[[matDimerName]] <- l[[matDimerName]] - transfer
    return(l)
  }))
}


amyloidAggregateGeneration_aggregateDissaggregationIntoMonomer <- function(neurons, maxDecline, stability, matMonomerName, matAggregateCountName, matAggregateSumName, aggregateMaxSize){
  return(mclapply(neurons, function(l){ #aMonomer has a single NA and aAggregateSum has a single NA
    if(sum(l[[matAggregateSumName]]) != 0){
      disAggMat <- rbinom(n = length(l[[matAggregateSumName]]),
                          size = l[[matAggregateSumName]],
                          prob = amyloidAggregateGeneration_logCurveMean(a = l[[matAggregateCountName]], b = l[[matAggregateSumName]], maxDecline = maxDecline, stability = stability, aggregateMaxSize = aggregateMaxSize))
      l[[matAggregateSumName]] <- l[[matAggregateSumName]] - disAggMat
      l[[matMonomerName]] <- l[[matMonomerName]] + disAggMat
    }
    return(l)
  }))
}


amyloidAggregateGeneration_logCurveMean <- function(a,b, maxDecline, stability, aggregateMaxSize){
  index <- a != 0
  m <- matrix(0, nrow(a), ncol(a))
  m[index] <- maxDecline - (maxDecline / (1 + exp(stability*(b[index]/aggregateMaxSize/a[index]) * (.5 + (b[index]/aggregateMaxSize/a[index])))))
  return(m)
}

amyloidAggregateGeneration_aggregateDeseeding <- function(neurons, maxDecline, stability, matMonomerName, matAggregateCountName, matAggregateSumName){
  #currently no de-seeding, until I find / remember reference
}
# 
# 
# logCurveMean <- function(a, b, maxDecline, stability){
#   mean(maxDecline / (1 + exp((stability / (a * 20 / b)) * (1 / (seq(b) / b) - (1 - (b / (a * 20)))))))
# }
# logCurveMean2 <- function(a,b, maxDecline, stability){
#   maxDecline / (1 + exp(stability*(b/20/a) * (.5 + (b/20/a))))
# }
# logCurveMean3 <- function(a,b, maxDecline, stability){
#   return(maxDecline / (1 + exp(stability*(b/20/a+1) * (.5 + (b/20/a+1)))))
# }
# 
# #L/(1 + exp(-k*(x - x_zero)))
# # 
# 
# a <- 1:20
# b <- 3:400
# dat <- expand.grid(a,b)
# dat$dep <- dat$Var1 / dat$Var2
# dat$ind <- 0
# dat$ind2 <- 0
# dat$ind3 <- 0
# for(i in 1:nrow(dat)){dat$ind[i] <- logCurveMean(dat$Var1[i], dat$Var2[i], maxDecline, stability)}
# for(i in 1:nrow(dat)){dat$ind2[i] <- logCurveMean2(dat$Var1[i], dat$Var2[i], maxDecline, stability)}
# for(i in 1:nrow(dat)){dat$ind3[i] <- logCurveMean3(dat$Var1[i], dat$Var2[i], maxDecline, stability)}
# plot(dat$ind2, dat$ind3)


# b <- 2000 #actual value does not matter
# maxDecline <- .3 # [0 - 1], real max Decline seems to be maxDecline/2 as the inflections point for the highest declining curve is at b
# stability <- 3 # (0 - Info) stability, hard to describe, but the range of 2-6 seems reasonable, at higher or lower values, the dependencies on a and b decline drastically
# 
# a <- floor(b/3)
# plot(logCurve(a,b, maxDecline, stability), col = "red")
# a <- floor(b/9)
# points(logCurve(a,b, maxDecline, stability), col ="blue")
# a <- floor(b/14)
# points(logCurve(a,b, maxDecline, stability), col ="green")
# a <- floor(b/20)
# points(logCurve(a,b, maxDecline, stability), col ="purple")
# curve3d(logCurveMean2(x,y,maxDecline,stability), from= c(1,1), to = c(20, 400), n = c(15, 150), sys3d = "image")
# curve3d(logCurveMean2(x,y,maxDecline,stability), from= c(1,1), to = c(20, 400), n = c(15, 150), sys3d = "contour", add =TRUE)

