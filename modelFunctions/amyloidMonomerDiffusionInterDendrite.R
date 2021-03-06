##########################################################
# This function diffuses amyloid monomers between 
# dendrite. The monomers are spread normally distributed up 
#and down a given number of dendrites. The function
# works in a way that it wraps arround the ends of the dendrite
# vector, s.t. its a ring vector
##########################################################

amyloidMonomerDiffusionInterDendrite <- function(neurons, factor, spread, spreadSD = .6, spreadMaxMultiplyer = 3, matName = "aMonomer"){
  require(parallel)
  parallel::mclapply(neurons, function(l){
    l[[matName]] <- amyloidMonomerDiffusionInterDendrite_VectorMultiplication(l[[matName]], factor, spread, spreadSD, spreadMaxMultiplyer)
    return(l)
  }) #, mc.preschedule = TRUE, mc.cores = 2
}

amyloidMonomerDiffusionInterDendrite_VectorMultiplication <- function(m, factor, spread, spreadSD, spreadMaxMultiplyer){
  shiftMat <- amyloidMonomerDiffusionInterDendrite_ShiftMatrix(dim(m)[1], factor, spread, spreadSD, spreadMaxMultiplyer)
  return(apply(m, 2, function(v) {
    nVec <- floor(shiftMat %*% v)
    d <- sum(v) - sum(nVec)
    i <- sample.int(n = length(v), size = d, replace = FALSE)
    nVec[i] <-  nVec[i] + 1
    return(nVec)
  }))
}

amyloidMonomerDiffusionInterDendrite_ShiftMatrix <- function(size, factor, spread, spreadSD, spreadMaxMultiplyer){ #maybe faster?
  if(factor >= 1) stop("spread / factor to large: You tried to spread all / more then available monomers.")
  if(spread  <= 0) stop("Spread of 0 (or negative) not possible!")
  if(spreadSD  < .1) stop("Standard deviation of < .1 not possible! This would introduce NaN (result to small).")
  if(size <= 1){stop("Matrix of size of 1 or smaller.")}
  if(spread > (size - 1) /2) {spread <- floor((size - 1) /2)}
  if(spread == 0){return(diag(1, size))}
  m <- diag(size)
  r <- runif(1, 1/spreadMaxMultiplyer, spreadMaxMultiplyer)
  sf <- amyloidMonomerDiffusionInterDendrite_spreadFactor(spread, factor * r, spreadSD)
  # sf <- sf * c(rep(r, spread), 1 - ((r - 1) * (1 / ((1 / factor) - 1))), rep(r, spread))
  for(i in 1:size){
    m[i, ((i-spread-1):(i+spread-1) %% size)+1] <- sf
  }
  return(m)
  
}

amyloidMonomerDiffusionInterDendrite_spreadFactor <- function(spread, factor, spreadSD = 1){
  v <- rep(1, 2 * spread + 1)
  i <- 1:spread
  i <- i / (spread / 2)
  s <- dnorm(i, sd = spreadSD)
  s <- s / sum(s) / 2 * factor
  v[1:spread] <- rev(s)
  v[(spread+2):(2*spread + 1)] <- s
  v[spread+1] <- 1 - factor
  if(any(is.na(v))){print(paste(spread, factor, spreadSD ,sep = " / "))}
  return(v)
}

# system.time(neurons <- amyloidMonomerDiffusionInterDendrite(neurons, .2, 1)) #takes 4 seconds
# system.time(test <- amyloidMonomerDiffusionInterDendrite_hack(neurons, .01, 5))
# source('http://www.stat.cmu.edu/~nmv/setup/mclapply.hack.R')

# ##### Example: #####
# mat <- matrix(0, nrow = 100, ncol = 100)
# mat[2,] <- 1:100
# mat
# (shiftMat <- amyloidMonomerDiffusionInterDendrite_ShiftMatrix(100, factor = .48, spread = 1, spreadSD = 1, spreadMaxMultiplyer = 2))
# library(plot.matrix)
# plot(amyloidMonomerDiffusionInterDendrite_VectorMultiplication(mat, factor = .3, spread = 1, spreadSD = 1, spreadMaxMultiplyer = 2))

