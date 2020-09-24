##########################################################
# This function diffuses amyloid monomers within the same 
# dendrite. The monomers are spread normally distributed to
# the left and to the right of the origin synapse. The function
# works in a way that it wraps arround the ends of the dendrite
# vector, s.t. its a ring vector
##########################################################

amyloidMonomerDiffusionIntraDendrite <- function(neurons, factor, spread, spreadSD = .6, spreadMaxMultiplyer = 3, matName = "aMonomer"){
  require(parallel)
  parallel::mclapply(neurons, function(l){
    l[[matName]] <- amyloidMonomerDiffusionIntraDendrite_VectorMultiplication(l[[matName]], factor, spread, spreadSD, spreadMaxMultiplyer)
    return(l)
  }) #, mc.preschedule = TRUE, mc.cores = 2
}

amyloidMonomerDiffusionIntraDendrite_VectorMultiplication <- function(m, factor, spread, spreadSD, spreadMaxMultiplyer){
  shiftMat <- amyloidMonomerDiffusionIntraDendrite_ShiftMatrix(dim(m)[2], factor, spread, spreadSD, spreadMaxMultiplyer)
  return(t(apply(m, 1, function(v) {
    nVec <- floor(shiftMat %*% v)
    d <- sum(v) - sum(nVec)
    i <- sample.int(n = length(v), size = d, replace = FALSE)
    nVec[i] <-  nVec[i] + 1
    return(nVec)
  })))
}

# amyloidMonomerDiffusionIntraDendrite_VectorMultiplicationSweep <- function(m, factor, spread, spreadSD, spreadMaxMultiplyer){ # SLOWER
#   shiftMat <- amyloidMonomerDiffusionIntraDendrite_ShiftMatrix(dim(m)[2], factor, spread, spreadSD, spreadMaxMultiplyer)
#   mat <- floor(sweep(m, 2, shiftMat, `%*%`, check.margin = FALSE))
#   d <- rowSums(m) - rowSums(mat)
#   colNumber <- ncol(m)
#   lapply(seq_along(d), function(i){
#     index <- sample.int(n = colNumber, size = d, replace = FALSE)
#     mat[i, index] <-  mat[i, index] + 1
#   })
#   return(mat)
# }

amyloidMonomerDiffusionIntraDendrite_ShiftMatrix <- function(size, factor, spread, spreadSD, spreadMaxMultiplyer){ #maybe faster?
  if(factor >= 1) stop("spread / factor to large: You tried to spread all / more then available monomers.")
  if(spread  == 0) stop("Spread of 0 not possible!")
  if(spreadSD  < .1) stop("Standard deviation of < .1 not possible! This would introduce NaN (result to small).")
  if(size <= 1) stop("Matrix of size of 1 or smaller.")
  if(spread > (size - 1) /2) {warning("Spread larger than Dendrite Size. Reducing spread size."); spread <- floor((size - 1) /2)}
  m <- matrix(0, nrow = size, ncol = size)
  r <- runif(1, 1/spreadMaxMultiplyer, spreadMaxMultiplyer)
  sf <- amyloidMonomerDiffusionIntraDendrite_spreadFactor(spread, min(factor * r, .999), spreadSD)
  # sf <- sf * c(rep(r, spread), 1 - ((r - 1) * (1 / ((1 / factor) - 1))), rep(r, spread))
  for(i in 1:size){
    m[i, ((i-spread-1):(i+spread-1) %% size)+1] <- sf
  }
  return(m)
}

amyloidMonomerDiffusionIntraDendrite_spreadFactor <- function(spread, factor, spreadSD = 1){
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

# system.time(neurons <- amyloidMonomerDiffusionIntraDendrite(neurons, .01, 90)) #takes 4 seconds
# system.time(test <- amyloidMonomerDiffusionIntraDendrite_hack(neurons, .01, 5))
# source('http://www.stat.cmu.edu/~nmv/setup/mclapply.hack.R')

# ##### Example: #####
# mat <- matrix(0, nrow = 100, ncol = 100)
# diag(mat) <- 1000
# mat[66,33] <- 1000
# mat
# # (shiftMat <- amyloidMonomerDiffusionIntraDendrite_ShiftMatrix(100, factor = .2, spread = 2, spreadSD = 1, spreadMaxMultiplyer = 2))
# # library(plot.matrix)
# # amyloidMonomerDiffusionIntraDendrite_VectorMultiplication(mat, factor = .5, spread = 50, spreadSD = 10, spreadMaxMultiplyer = 2)
# library(reshape2)
# library(ggplot2)
# matSpreaded <- amyloidMonomerDiffusionIntraDendrite_VectorMultiplication(mat, factor = .9, spread = 20, spreadSD = .6, spreadMaxMultiplyer = 1)
# conMatLong <- melt(matSpreaded)
# ggplot(conMatLong) +
#   geom_raster(aes(y = abs(Var1 - 101), x = Var2, fill = value))

