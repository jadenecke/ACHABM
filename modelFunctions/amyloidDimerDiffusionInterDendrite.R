##########################################################
# This function diffuses amyloid Dimers between 
# dendrite. The Dimers are spread normally distributed up 
# and down a given number of dendrites. The function
# works in a way that it wraps arround the ends of the dendrite
# vector, s.t. its a ring vector
##########################################################

amyloidDimerDiffusionInterDendrite <- function(neurons, factor, spread, spreadSD = .6, spreadMaxMultiplyer = 3, matName = "aDimer"){
  require(parallel)
  parallel::mclapply(neurons, function(l){
    l[[matName]] <- amyloidDimerDiffusionInterDendrite_VectorMultiplication(l[[matName]], factor, spread, spreadSD, spreadMaxMultiplyer)
    return(l)
  }) #, mc.preschedule = TRUE, mc.cores = 2
}

amyloidDimerDiffusionInterDendrite_VectorMultiplication <- function(m, factor, spread, spreadSD, spreadMaxMultiplyer){
  shiftMat <- amyloidDimerDiffusionInterDendrite_ShiftMatrix(dim(m)[1], factor, spread, spreadSD, spreadMaxMultiplyer)
  return(apply(m, 2, function(v) {
    nVec <- floor(shiftMat %*% v)
    d <- sum(v) - sum(nVec)
    i <- sample.int(n = length(v), size = d, replace = FALSE)
    nVec[i] <-  nVec[i] + 1
    return(nVec)
  }))
}

amyloidDimerDiffusionInterDendrite_ShiftMatrix <- function(size, factor, spread, spreadSD, spreadMaxMultiplyer){ #maybe faster?
  if(factor >= 1) stop("spread / factor to large: You tried to spread all / more then available Dimers.")
  if(spread  <= 0) stop("Spread of 0 (or negative) not possible!")
  if(spreadSD  < .1) stop("Standard deviation of < .1 not possible! This would introduce NaN (result to small).")
  if(size <= 1){stop("Matrix of size of 1 or smaller.")}
  if(spread > (size - 1) /2) {warning("Spread larger than number of Dendrites. Reducing spread size."); spread <- floor((size - 1) /2)}
  m <- diag(size)
  r <- runif(1, 1/spreadMaxMultiplyer, spreadMaxMultiplyer)
  sf <- amyloidDimerDiffusionInterDendrite_spreadFactor(spread, factor * r, spreadSD)
  # sf <- sf * c(rep(r, spread), 1 - ((r - 1) * (1 / ((1 / factor) - 1))), rep(r, spread))
  for(i in 1:size){
    m[i, ((i-spread-1):(i+spread-1) %% size)+1] <- sf
  }
  return(m)
  
}

amyloidDimerDiffusionInterDendrite_spreadFactor <- function(spread, factor, spreadSD = 1){
  v <- rep(1, 2 * spread + 1)
  i <- 1:spread
  i <- i / (spread / 2)
  s <- dnorm(i, sd = spreadSD)
  s <- s / sum(s) / 2 * factor
  v[1:spread] <- rev(s)
  v[(spread+2):(2*spread + 1)] <- s
  v[spread+1] <- 1 - factor
  return(v)
}

# system.time(neurons <- amyloidDimerDiffusionInterDendrite(neurons, .2, 1)) #takes 4 seconds
# system.time(test <- amyloidDimerDiffusionInterDendrite_hack(neurons, .01, 5))
# source('http://www.stat.cmu.edu/~nmv/setup/mclapply.hack.R')

###### Example: #####
# mat <- matrix(1, nrow = 3, ncol = 10)
# mat[2,] <- 1:10
# mat
# (shiftMat <- amyloidDimerDiffusionInterDendrite_ShiftMatrix(4, factor = .48, spread = 1, spreadSD = 1, spreadMaxMultiplyer = 2))
# amyloidDimerDiffusionInterDendrite_VectorMultiplication(mat, factor = .3, spread = 1, spreadSD = 1, spreadMaxMultiplyer = 2)

