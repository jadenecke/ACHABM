##########################################################
# This function takes a neuron list and 
# updates the amyloid monomer matrices of every neuron
##########################################################

amyloidMonomerGeneration <- function(neurons, mean, sd, upperBound, matName = "aMonomer", activityName = "activity"){
  #if(all(unlist(lapply(neurons, function(l) (is.null(l[[matName]]) || is.null(l[[activityName]])))))){stop(paste("Can't apply amyloidMonomerGeneration, either", matName, "or", activityName, "missing!"))}
  for(i in seq_along(neurons)){
    if(neurons[[i]][["alive"]]){
      neurons[[i]][[matName]] <- neurons[[i]][[matName]] + pmax(pmin(round(rnorm(n = length(neurons[[i]][[matName]]), mean, sd) * neurons[[i]][[activityName]]), upperBound), 0)
    }
  }
  return(neurons)
}


# updatedAmyloidMonomers2 <- function(neurons, d = runif, matName = "aMonomer", ... ){
#   lapply(seq_along(neurons), function(i){
#     neurons[[i]][[matName]] <- neurons[[i]][[matName]] + d(length(neurons[[i]][[matName]]), ...)
#     neurons[[i]]
#   })
# }

#  # microbenchmark::microbenchmark(list.update(neurons, z = z + 1), lapply(seq_along(neurons), function(i){neurons[[i]][["z"]] <- neurons[[i]][["z"]] + 1; neurons[[i]]})) #list.update function is faster
# uz <- function(neurons){ # for-loop function fastest!
#   for(i in seq_along(neurons)){
#     neurons[[i]]$z <<- neurons[[i]]$z + 1 ## <<- changes global variable: CAREFULL! this might change the result of subsequent calculations!
#   }
#   return(neurons)
# }