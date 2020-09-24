##########################################################
# spreads amyloid Dimers from one synapse to its corresponding 
# neurons, where it is then distributed equally over all synasps.
# The probability to spread Dimer is dependent on the amount
# of Dimers available in a logistic fashion
##########################################################

amyloidDimerDiffusionInterNeuron <- function(neurons, matName = "aDimer", stackName = "aDimerStack", connectionsName = "conMat", maximumSpreadProbability, spreadDependencyCurveSteepness, spreadDependencyCurveInflectionPoint){
  require(parallel)
  spreadList <- amyloidDimerDiffusionInterNeuron_calculateAmountToTransfer(neurons, matName, maximumSpreadProbability, spreadDependencyCurveSteepness, spreadDependencyCurveInflectionPoint)
  neurons <- amyloidDimerDiffusionInterNeuron_transferFromNeuronToNeuron(neurons, spreadList, matName, stackName, connectionsName)
  neurons <- amyloidDimerDiffusionInterNeuron_spreadFromNeuronToSynapses(neurons, matName, stackName)
  return(neurons)
}

amyloidDimerDiffusionInterNeuron_calculateAmountToTransfer <- function(neurons, matName, maximumSpreadProbability, spreadDependencyCurveSteepness, spreadDependencyCurveInflectionPoint){
  return(parallel::mclapply(neurons, function(l){
    if(l[["alive"]]){
    spreadVec <- rbinom(length(l[[matName]]),
                        l[[matName]],
                        maximumSpreadProbability/(1 + exp(-spreadDependencyCurveSteepness*(l[[matName]] - spreadDependencyCurveInflectionPoint))))
    return(spreadVec)
    } else {
      return(c(0, 0))
    }
  })) #, mc.preschedule = TRUE, mc.cores = 2
}
# curve(L/(1 + exp(-k*(x - x_zero))), 1, 10000) # L = Maxima; k = steepness; x_zero = inflection Point (Wendepunkt)
# curve(maximumSpreadProbability/(1 + exp(-spreadDependencyCurveSteepness*(x - spreadDependencyCurveInflectionPoint))), 1, 1000) 

amyloidDimerDiffusionInterNeuron_transferFromNeuronToNeuron <- function(neurons, spreadList, matName, stackName, connectionsName){
  for(i in seq(length(spreadList))){
    if(sum(spreadList[[i]]) != 0){
      neurons[[i]][[matName]] <- neurons[[i]][[matName]] - spreadList[[i]]
      con <- neurons[[i]][[connectionsName]]
      for(j in seq_along(con)){
        neurons[[con[j]]][[stackName]] <- neurons[[con[j]]][[stackName]] + spreadList[[i]][j]
      }
    }
  }
  return(neurons)
}

amyloidDimerDiffusionInterNeuron_spreadFromNeuronToSynapses <- function(neurons, matName, stackName){
  return(lapply(neurons, function(l){
    if(l[[stackName]] != 0 & l[["alive"]]){
      s <- sample(amyloidDimerDiffusionInterNeuron_spreadVec(l[[stackName]], length(l[[matName]])))
      l[[matName]] <- l[[matName]] + s
      l[[stackName]] <- 0
    }
    return(l)
  })) #, mc.preschedule = TRUE, mc.cores = 2
}

amyloidDimerDiffusionInterNeuron_spreadVec <- function(n, l){
  x <- n / l
  whole = floor(x)
  fraction = x - whole
  return(rep(whole, l) + c(rep(1, round(l * fraction)), rep(0, l - round(l * fraction))))
}



