##########################################################
# Amyloid Dimer Disaggregaton function: follows a binomial
# distribution with a fixed probability. One dimer disaggregates
# into two monomers
##########################################################

amyloidDimerDisaggregation <- function(neurons, disaggregationProbability, matDimerName = "aDimer", matMonomerName = "aMonomer"){
  require("parallel")
  return(mclapply(neurons, function(l){
    genMat <- rbinom(n = length(l[[matDimerName]]),
                     size = l[[matDimerName]],
                     prob = disaggregationProbability
    )
    l[[matDimerName]] <- l[[matDimerName]]  - genMat
    l[[matMonomerName]] <- l[[matMonomerName]]  + (genMat * 2)
    return(l)
  }))
}
