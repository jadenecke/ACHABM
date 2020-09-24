##########################################################
# 
##########################################################

amyloidPlaqueGenerationFromAll <- function(neurons,
                                    plaqueMaximumSize,
                                    plaquePullIntraDendrite,
                                    plaquePullInterDendrite,
                                    plaquePullMaxProb,
                                    softLimit,
                                    AggregatesizePullRelation, #(0 , 1) values close to 0 indicate only large Aggregates are used to fill Plaques
                                    amyloidAggregateGeneration_aggregateMaxSize,
                                    plaqueSeedProbabilityCurveMax,
                                    plaqueSeedProbabilityCurveSteepness,
                                    plaqueSeedProbabilityInflectionPoint,
                                    aPlaqueName = "aPlaque",
                                    aAggregateSumName = "aAggregateSum",
                                    aAggregateCountName = "aAggregateCount",
                                    aDimerName = "aDimer"
){
  require("parallel")
  require("Matrix")
  return(mclapply(neurons, function(l){
    #grow existing plaques first (so that no new plaques form as long as there are plaques that can sequester aggregates):
    plaqueIndex <- which(0 < l[[aPlaqueName]] & l[[aPlaqueName]] < (plaqueMaximumSize - softLimit), arr.ind = TRUE)
    if(length(plaqueIndex) > 0){
      for(i in sample(seq(dim(plaqueIndex)[1]))){
        pull <- amyloidPlaqueGeneration_getPullMatrix(pos = plaqueIndex[i,],
                                                      d = dim(l[[aPlaqueName]]),
                                                      plaquePullIntraDendrite = plaquePullIntraDendrite, 
                                                      plaquePullInterDendrite = plaquePullInterDendrite
        )
        
        transferAggregate <- rbinom(sum(pull[["pullIndex"]]),
                           l[[aAggregateSumName]][pull[["pullIndex"]]],
                           pull[["pullFactor"]][pull[["pullIndex"]]] * plaquePullMaxProb
        )
        
        transferDimer <- rbinom(sum(pull[["pullIndex"]]),
                                    l[[aDimerName]][pull[["pullIndex"]]],
                                    pull[["pullFactor"]][pull[["pullIndex"]]] * plaquePullMaxProb
        )
        
        transferAggregate[transferAggregate < 3] <- 0
        maxTransfer <- plaqueMaximumSize - l[[aPlaqueName]][plaqueIndex[i,1], plaqueIndex[i,2]]

        if((sum(transferAggregate) - maxTransfer)  > 0){
          transferAggregate <- transferAggregate - ceiling((sum(transferAggregate) - maxTransfer ) / length(transferAggregate))
        }
        
        l[[aAggregateSumName]][pull[["pullIndex"]]] <- l[[aAggregateSumName]][pull[["pullIndex"]]] - transferAggregate
        l[[aPlaqueName]][plaqueIndex[i,1], plaqueIndex[i,2]] <- l[[aPlaqueName]][plaqueIndex[i,1], plaqueIndex[i,2]] + sum(transferAggregate) 
        
        aMax <-  pmin(floor(transferAggregate / 3), floor(l[[aAggregateCountName]][pull[["pullIndex"]]] * plaquePullMaxProb))
        aMin <- ceiling(transferAggregate / amyloidAggregateGeneration_aggregateMaxSize)
        maxReduce <- ceiling(l[[aAggregateSumName]][pull[["pullIndex"]]] / amyloidAggregateGeneration_aggregateMaxSize)
        l[[aAggregateCountName]][pull[["pullIndex"]]] <- pmax(maxReduce, l[[aAggregateCountName]][pull[["pullIndex"]]] - (rbinom(length(transferAggregate),
                                                                                                                                 pmax(0, aMax - aMin),
                                                                                                                                 AggregatesizePullRelation)
                                                                                                                          + aMin)
        )
        
        l[[aDimerName]][pull[["pullIndex"]]] <- l[[aDimerName]][pull[["pullIndex"]]] - transferDimer
        l[[aPlaqueName]][plaqueIndex[i,1], plaqueIndex[i,2]] <- l[[aPlaqueName]][plaqueIndex[i,1], plaqueIndex[i,2]] + sum(transferDimer) * 2 
      }
    }
    
    #seed new plaques:

    if(any(l[[aAggregateCountName]] != 0)){
      if(ceiling(plaquePullIntraDendrite / 3) >= (dim(l[[aPlaqueName]])[2] - 1)){
        spreadMat <- matrix(1, nrow = dim(l[[aPlaqueName]])[2], ncol = dim(l[[aPlaqueName]])[2])
      } else {
        spreadMat <- band(matrix(1, dim(l[[aPlaqueName]])[2], dim(l[[aPlaqueName]])[2]), -ceiling(plaquePullIntraDendrite / 3), ceiling(plaquePullIntraDendrite / 3))
        
      }
      
      l[[aPlaqueName]] <- l[[aPlaqueName]] + rbinom(n = length(l[[aPlaqueName]]),
                                                    size = as.matrix(l[[aPlaqueName]] %*% spreadMat) == 0,
                                                    prob = plaqueSeedProbabilityCurveMax/(1 + exp(-plaqueSeedProbabilityCurveSteepness * ((l[[aAggregateCountName]] +  l[[aDimerName]]) - plaqueSeedProbabilityInflectionPoint)))
      )
    }
    return(l)
  }))
}

# curve(1/(1 + exp(-1 * (.001 - (1/x)))), 1,50, ylim = c(0,1))



# m <- matrix(0, ncol = 2000, 2000)
# m[getPullIndex(22217255, nrow(m), ncol(m), 500, 500)] <- 1
# plot(raster(m))

# d <- data.frame(i = floor((1:21)/2))
# d$sin <- sin(pi/2*(1-(floor((1:21)/2))/(10 + 1)))
# d$cos <- cos(pi/2*(1-(floor((1:21)/2))/(10 + 1)))
# d$diff <- i - d$cos
# View(d)
# 
# plot(sin(pi/2*(1-(floor((
#   floor((1:21)/2) * cos(pi/2*(1-(floor((1:21)/2))/(10 + 1)))
#   )/2))/(10 + 1))))
# 
# curvePoints <- function(n){
#   n_full <- (n + 1) *4
#   p <- seq(n_full)
#   x <- tan(seq(1, -1 , length.out = n_full))
#   a <- (p*x)  * (360 / (n_full))
#   s <- sin(deg2rad(a[seq(n)]))
#   return(sort(c(1, s, s), decreasing = TRUE))
# }
# 
# curvePoints <- function(n){
#   p <- seq(n*4)
#   a <- p * (360 / (n*4))
#   return(data.frame(x = cos(deg2rad(a[1:3])),
#                     y = sin(deg2rad(a[1:3]))))
# }
# 
# View(curvePoints(3))
# plot(curvePoints(10))
# 
# rad2deg <- function(rad) {(rad * 180) / (pi)}
# deg2rad <- function(deg) {(deg * pi) / (180)}

# N=400; M=400;
# cx=200; cy=100
# radiusX=50; radiusY = 50
# 
# binaryMap=outer((x-cx)^2, ((y*2)-cy)^2, "+") <= radius^2
# plot(raster(binaryMap))

amyloidPlaqueGeneration_getPullMatrix <- function(pos, d, plaquePullIntraDendrite, plaquePullInterDendrite){
  if(plaquePullInterDendrite > plaquePullIntraDendrite){
    f <- plaquePullInterDendrite / plaquePullIntraDendrite
    x <- amyloidPlaqueGeneration_circShift(amyloidPlaqueGeneration_getVec(d[1]), ceiling(d[1] / 2 - pos[1]))
    y <- amyloidPlaqueGeneration_circShift(amyloidPlaqueGeneration_getVec(d[2]), ceiling(d[2] / 2 - pos[2]))
    l <- list(pullIndex = outer((x)^2, ((y)*f)^2, "+") <= plaquePullInterDendrite^2,
              pullFactor = matrix(pmax(0, (plaquePullInterDendrite^2 + 1 - outer((x)^2, ((y)*f)^2, "+"))/(plaquePullInterDendrite^2 + 1)), d[1], d[2])
    )
    return(l)
  } else if(plaquePullInterDendrite < plaquePullIntraDendrite){
    f <- plaquePullIntraDendrite / plaquePullInterDendrite
    x <- amyloidPlaqueGeneration_circShift(amyloidPlaqueGeneration_getVec(d[1]), ceiling(d[1] / 2 - pos[1]))
    y <- amyloidPlaqueGeneration_circShift(amyloidPlaqueGeneration_getVec(d[2]), ceiling(d[2] / 2 - pos[2]))
    l <- list(pullIndex = outer(((x) * f)^2, (y)^2, "+") <= plaquePullIntraDendrite^2,
              pullFactor = pmax(0, (plaquePullIntraDendrite^2 + 1 - outer(((x) * f)^2, (y)^2, "+"))/(plaquePullIntraDendrite^2 + 1))
    )
    return(l)
  } else if(plaquePullInterDendrite == plaquePullIntraDendrite){
    x <- amyloidPlaqueGeneration_circShift(amyloidPlaqueGeneration_getVec(d[1]), ceiling(d[1] / 2 - pos[1]))
    y <- amyloidPlaqueGeneration_circShift(amyloidPlaqueGeneration_getVec(d[2]), ceiling(d[2] / 2 - pos[2]))
    l <- list(pullIndex = outer((x)^2, (y)^2, "+") <= plaquePullInterDendrite^2,
              pullFactor = pmax(0, (plaquePullInterDendrite^2 + 1 - outer((x)^2, (y)^2, "+"))/(plaquePullInterDendrite^2 + 1))
    )
    return(l)
  } else {
    stop("impossible")
    return(NULL)
  }
}

amyloidPlaqueGeneration_getVec <- function(len){
  if(len %% 2 == 1){
    n <- floor(len/2)
    return(c(seq(n, 1), 0 , seq(n)))
  } else {
    n <- floor(len/2)
    return(c(seq(n-1, 1), 0 , seq(n)))
  }
}

amyloidPlaqueGeneration_circShiftSlower <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

amyloidPlaqueGeneration_circShift <- function(x, n = 1) {
  if (n == 0) {
    return(x)
  } else if(n > 0) {
    return(c(x[(n+1):length(x)], x[1:n]))
  } else if(n < 0){
    l <- length(x)
    return(c(x[(l - abs(n) + 1):l], x[1:(l - abs(n))]))
  } else {
    stop("something fucked up with amyloidPlaqueGeneration_circShift")
    return(NULL)
  }
}

# pm <- matrix(0, 200, 200)
# pm[20,80] <- 1
# pm[130,130] <- 1
# pm
# a <- matrix(100, 200, 200)
# plaqueIndex <- which(pm != 0, arr.ind = TRUE)
# for(i in sample(seq(dim(plaqueIndex)[1]))){
#   pull <- amyloidPlaqueGeneration_getPullMatrix(pos = plaqueIndex[i,],
#                                                 d = dim(pm),
#                                                 plaquePullIntraDendrite = 70, 
#                                                 plaquePullInterDendrite = 50
#   )
#   transferAggregate <- rbinom(sum(pull[["pullIndex"]]),
#                      a[pull[["pullIndex"]]],
#                      pull[["pullFactor"]][pull[["pullIndex"]]] * .8
#   )
#   a[pull[["pullIndex"]]] <- a[pull[["pullIndex"]]] - transferAggregate
#   pm[plaqueIndex[i,1], plaqueIndex[i,2]] <- pm[plaqueIndex[i,1], plaqueIndex[i,2]] + sum(transferAggregate)
# }
# 
# #plot(pm)
# plot(a)
# 
# 
# m <- matrix(0, 6,10)
# m[3,5] <- 1
# mt <- diag(1, nrow = dim(m)[2], ncol = dim(m)[2])
# diag(mt[-1,]) <- 1
# diag(mt[,-1]) <- 1
# m %*% mt

