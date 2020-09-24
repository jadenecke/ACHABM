##########################################################
# Amyloid dimer generation: Two monomers are transformed into
# one dimer with the amount of monomers transfered determined
# by a binomial distribution with its probability given by logistic distribution
##########################################################

amyloidDimerGeneration <- function(neurons, matDimerName = "aDimer", matMonomerName = "aMonomer", maximumPercentTransform, probCurveSteepness, probCurveInflectionPoint){
  require("parallel")
  return(mclapply(neurons, function(l){
    if(l[["alive"]]){
    genMat <- rbinom(length(l[[matMonomerName]]),
                     l[[matMonomerName]],
                     maximumPercentTransform/(1 + exp(-probCurveSteepness*(l[[matMonomerName]] - probCurveInflectionPoint)))
                     )
    genMat <- floor((genMat) / 2)
    l[[matDimerName]] <- l[[matDimerName]]  + genMat
    l[[matMonomerName]] <- l[[matMonomerName]]  - (genMat * 2)
    }
    return(l)
  }))
}
