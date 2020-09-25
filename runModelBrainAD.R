rm(list = ls())

# Dependencies
library(rlist)
library(parallel)
source("modelParameter.R")
source("initiateModel.R")
source("generalFunctions.R")
lapply(list.files("modelFunctions",full.names = TRUE), source)


#make Paths and logging stuff:
logPath <- file.path("G:", "Uni", "SEDS MA Data", "logs")
logPath <- file.path(logPath, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_brainAD"))
dir.create(logPath)
logParamter(logPath =  logPath)

#generate dummy network:
# source(file.path("networkGeneration", "genN-ClusterNetwork.R"))
# neurons <- genNClusterNetwork(numberOfCluster = numberOfCluster,
#                               clusterMinSize = clusterMinSize,
#                               clusterMaxSize = clusterMaxSize,
#                               dendritesPerNeuronMin = dendritesPerNeuronMin,
#                               dendritesPerNeuronMax = dendritesPerNeuronMax,
#                               synapsesPerDendriteMin = synapsesPerDendriteMin,
#                               synapsesPerDendriteMax = synapsesPerDendriteMax)
# 

#generate Beaujoin 2018 Model
source(file.path("networkGeneration", "genBeaujoin2018.R"))
neurons <- initiateModel(size = size,
                         connectionMatrix = conMat,
                         dendritesPerNeuronMin = 1,
                         dendritesPerNeuronMax = 4,
                         synapsesPerDendriteMin = 7,
                         synapsesPerDendriteMax = 70)

# hist(unlist(lapply(neurons, function(l){nrow(l$conMat)})))
# table(unlist(lapply(neurons, function(l){ncol(l$conMat)})))
# conMat %*% diag(1/ size)


#set up logging results:
sumAmyloidMonomers <- rep(0, numberOfSteps)
sumAmyloidDimers <- rep(0, numberOfSteps)
sumAmyloidAggregatesCount <- rep(0, numberOfSteps)
sumAmyloidAggregatesSize <- rep(0, numberOfSteps)
sumAmyoidPlaques <- rep(0, numberOfSteps)
countAmyoidPlaques <- rep(0, numberOfSteps)
synapseLoss <- rep(0, numberOfSteps)
synapseActivity <- rep(0, numberOfSteps)
sumNFT <- rep(0, numberOfSteps)
meanNFTSeedProbability <- rep(0, numberOfSteps)
meanNeuronActivity <- rep(0, numberOfSteps)
alivePerc <- rep(0, numberOfSteps)
numberOfSynapses <- sum(unlist(lapply(neurons, function(l) length(l$aMonomer))))


# #burn in:
# lastStepMonomers <- 1
# currentStepMonomers <- 2
# i <- 0
# burnInLogList <- list()
# while(currentStepMonomers / lastStepMonomers > 1.1){
#   i <- i + 1
#   print(paste0("-------- Burn in Step ", i, " --------"))
#   print(paste0("Current Ratio: ", currentStepMonomers / lastStepMonomers ))
#   print(paste0("Generating Amyloid Monomers"))
#   neurons <- amyloidMonomerGeneration(neurons,
#                                       mean = amyloidMonomerGeneration_mu,
#                                       sd = amyloidMonomerGeneration_delta,
#                                       upperBound = amyloidMonomerGeneration_upperBound
#   )
#   print(paste0("Cleaning Amyloid"))
#   neurons <- amyloidMonomerCleaning(neurons,
#                                     mean = amyloidMonomerCleaning_mu,
#                                     sd = amyloidMonomerCleaning_delta,
#                                     upperBound = amyloidMonomerCleaning_upperBound,
#                                     currentStep = 1,
#                                     declineFactor = amyloidMonomerCleaning_declineFactor,
#                                     declineFlag = FALSE
#   )
#   lastStepMonomers <- currentStepMonomers
#   currentStepMonomers <- sum(unlist(lapply(neurons, function(x) x$aMonomer))) / numberOfSynapses
#   burnInLogList <- list.append(burnInLogList, currentStepMonomers)
# 
# }
# plot(unlist(burnInLogList))


#run Model:
for(i in 1:numberOfSteps){
  print(paste0("-------- Step ", i, " / ", numberOfSteps, " --------"))
  print(paste0("Generating Amyloid Monomers"))
  neurons <- amyloidMonomerGeneration(neurons,
                                      mean = amyloidMonomerGeneration_mu,
                                      sd = amyloidMonomerGeneration_delta,
                                      upperBound = amyloidMonomerGeneration_upperBound,
                                      matName = "aMonomer",
                                      activityName = "activity"
  )
  print(paste0("Cleaning Amyloid"))
  neurons <- amyloidMonomerCleaning(neurons,
                                    mean = amyloidMonomerCleaning_mu,
                                    sd = amyloidMonomerCleaning_delta,
                                    upperBound = amyloidMonomerCleaning_upperBound,
                                    currentStep = i,
                                    declineFactor = amyloidMonomerCleaning_declineFactor,
                                    declineFlag = TRUE,
                                    matName = "aMonomer",
                                    activityName = "activity"
  )
  
  print(paste0("Diffusing amyloid Monomers - Intra Dendrite Diffusion"))
  neurons <- amyloidMonomerDiffusionIntraDendrite(neurons = neurons,
                                                  spread = amyloidMonomerDiffusionIntraDendrite_range,
                                                  factor = amyloidMonomerDiffusionIntraDendrite_amountOfDiffusion,
                                                  spreadSD = amyloidMonomerDiffusionIntraDendrite_spreadSD,
                                                  spreadMaxMultiplyer = amyloidMonomerDiffusionIntraDendrite_spreadMaxMultiplyer,
                                                  matName = "aMonomer"
  )
  
  print(paste0("Diffusing amyloid Monomers - Inter Dendrite Diffusion"))
  neurons <- amyloidMonomerDiffusionInterDendrite(neurons = neurons,
                                                  spread = amyloidMonomerDiffusionInterDendrite_range,
                                                  factor = amyloidMonomerDiffusionInterDendrite_amountOfDiffusion,
                                                  spreadSD = amyloidMonomerDiffusionInterDendrite_spreadSD,
                                                  spreadMaxMultiplyer = amyloidMonomerDiffusionInterDendrite_spreadMaxMultiplyer,
                                                  matName = "aMonomer"
  )
  
  print(paste0("Diffusing amyloid Monomers - Inter Neuron Diffusion"))
  neurons <- amyloidMonomerDiffusionInterNeuron(neurons = neurons,
                                                matName = "aMonomer",
                                                stackName = "aMonomerStack",
                                                connectionsName = "conMat",
                                                maximumSpreadProbability = amyloidMonomerDiffusionInterNeuron_maximumSpreadProbability,
                                                spreadDependencyCurveSteepness = amyloidMonomerDiffusionInterNeuron_spreadDependencyCurveSteepness,
                                                spreadDependencyCurveInflectionPoint = amyloidMonomerDiffusionInterNeuron_spreadDependencyCurveInflectionPoint
  )
  
  print(paste0("Generating Amyloid Dimers"))
  neurons <- amyloidDimerGeneration(neurons = neurons,
                                    matDimerName = "aDimer",
                                    matMonomerName = "aMonomer",
                                    maximumPercentTransform = amyloidDimerGeneration_maximumPercentTransform,
                                    probCurveSteepness = amyloidDimerGeneration_probCurveSteepness,
                                    probCurveInflectionPoint = amyloidDimerGeneration_probCurveInflectionPoint
  )
  
  print(paste0("Diffusing amyloid Dimers - Intra Dendrite Diffusion"))
  neurons <- amyloidDimerDiffusionIntraDendrite(neurons = neurons,
                                                spread = amyloidDimerDiffusionIntraDendrite_range,
                                                factor = amyloidDimerDiffusionIntraDendrite_amountOfDiffusion,
                                                spreadSD = amyloidDimerDiffusionIntraDendrite_spreadSD,
                                                spreadMaxMultiplyer = amyloidDimerDiffusionIntraDendrite_spreadMaxMultiplyer,
                                                matName = "aDimer"
  )
  
  print(paste0("Diffusing amyloid Dimers - Inter Dendrite Diffusion"))
  neurons <- amyloidDimerDiffusionInterDendrite(neurons = neurons,
                                                spread = amyloidDimerDiffusionInterDendrite_range,
                                                factor = amyloidDimerDiffusionInterDendrite_amountOfDiffusion,
                                                spreadSD = amyloidDimerDiffusionInterDendrite_spreadSD,
                                                spreadMaxMultiplyer = amyloidDimerDiffusionInterDendrite_spreadMaxMultiplyer,
                                                matName = "aDimer"
  )
  
  print(paste0("Diffusing amyloid Dimers - Inter Neuron Diffusion"))
  neurons <- amyloidDimerDiffusionInterNeuron(neurons = neurons,
                                              matName = "aDimer",
                                              stackName = "aDimerStack",
                                              connectionsName = "conMat",
                                              maximumSpreadProbability = amyloidDimerDiffusionInterNeuron_maximumSpreadProbability,
                                              spreadDependencyCurveSteepness = amyloidDimerDiffusionInterNeuron_spreadDependencyCurveSteepness,
                                              spreadDependencyCurveInflectionPoint = amyloidDimerDiffusionInterNeuron_spreadDependencyCurveInflectionPoint
  )
  
  print(paste0("Amyloid Dimer Dissagregation"))
  neurons <- amyloidDimerDisaggregation(neurons = neurons,
                                        disaggregationProbability = amyloidDimerDisaggregation_DisaggregationProbability,
                                        matDimerName = "aDimer",
                                        matMonomerName = "aMonomer"
  )
  
  print(paste0("Amyloid Aggregate Generation"))
  neurons <- amyloidAggregateGeneration(neurons,
                                        seedingProbabiltyMax = amyloidAggregateGeneration_seedingProbabiltyMax,
                                        seedingProbabilityCurveSteepness = amyloidAggregateGeneration_seedingProbabilityCurveSteepness,
                                        seedingProbabilityCurveInflectionPoint = amyloidAggregateGeneration_seedingProbabilityCurveInflectionPoint,
                                        aggregateGrowthProbabiltyMax = amyloidAggregateGeneration_aggregateGrowthProbabiltyMax,
                                        aggregateGrowthProbabilityCurveSteepness = amyloidAggregateGeneration_aggregateGrowthProbabilityCurveSteepness,
                                        aggregateGrowthDelay = amyloidAggregateGeneration_aggregateGrowthDelay,
                                        maxDecline = amyloidAggregateGeneration_maxDecline,
                                        stability = amyloidAggregateGeneration_stability,
                                        aggregateMaxSize = amyloidAggregateGeneration_aggregateMaxSize,
                                        matMonomerName = "aMonomer",
                                        matDimerName = "aDimer",
                                        matAggregateCountName = "aAggregateCount",
                                        matAggregateSumName = "aAggregateSum"
  )
  
  print(paste0("Amyloid Plaque Generation & Growth"))
  neurons <- amyloidPlaqueGenerationFromAll(neurons, #reduced plaque max growth /3 because now we have 3 sources -> .005 to .0015
                                     plaqueMaximumSize = amyloidPlaqueGeneration_plaqueMaximumSize,
                                     plaquePullIntraDendrite = amyloidPlaqueGeneration_plaquePullIntraDendrite,
                                     plaquePullInterDendrite = amyloidPlaqueGeneration_plaquePullInterDendrite,
                                     plaquePullMaxProb = amyloidPlaqueGeneration_plaquePullMaxProb,
                                     AggregatesizePullRelation = amyloidPlaqueGeneration_AggregatesizePullRelation, #(0 , 1) values close to 0 indicate only large Aggregates are used to fill Plaques
                                     amyloidAggregateGeneration_aggregateMaxSize = amyloidAggregateGeneration_aggregateMaxSize,
                                     plaqueSeedProbabilityCurveMax = amyloidPlaqueGeneration_plaqueSeedProbabilityCurveMax,
                                     plaqueSeedProbabilityCurveSteepness = amyloidPlaqueGeneration_plaqueSeedProbabilityCurveSteepness,
                                     plaqueSeedProbabilityInflectionPoint = amyloidPlaqueGeneration_plaqueSeedProbabilityInflectionPoint,
                                     softLimit = amyloidPlaqueGeneration_softLimit,
                                     aPlaqueName = "aPlaque",
                                     aAggregateSumName = "aAggregateSum",
                                     aAggregateCountName = "aAggregateCount"
  )
  
  print(paste0("NFT Generation & Clearance"))
  neurons <- nftGeneration(neurons,
                           nftAcceleration = nftGeneration_nftAcceleration,
                           probabilitySaclingFactor = nftGeneration_probabilitySaclingFactor,
                           maxNftGrowth = nftGeneration_maxNftGrowth, 
                           synapseActivityUpdate_nftCutOffMean = synapseActivityUpdate_nftCutOffMean,
                           nftFlatClearance = nftGeneration_nftFlatClearance,
                           nftName = "nft",
                           aAggregateCountName = "aAggregateCount",
                           nftSeedProbability = "nftSeedProbability",
                           activityName = "activity"
  )
  
  print(paste0("NFT Seed Probability Spreading"))
  neurons <- nftSeedProbabilitySpread(neurons,
                                      lag = nftSeedProbabilitySpread_lag,
                                      maxSeedProb = nftSeedProbabilitySpread_maxSeedProb,
                                      seedProbMaxRise = nftSeedProbabilitySpread_seedProbMaxRise,
                                      nftName = "nft",
                                      nftSeedProbabilityName = "nftSeedProbability",
                                      nftSeedProbabilityStackName = "nftSeedProbabilityStack",
                                      connectionsName = "conMat"
  )
  
  print(paste0("Update Synapse Activity"))
  neurons <- synpaseActivityUpdate(neurons,
                                   declineStability = synapseActivityUpdate_declineStability,
                                   declineStart = synapseActivityUpdate_declineStart,
                                   activityCutOff = synapseActivityUpdate_activityCutOff,
                                   nftCutOffMean = synapseActivityUpdate_nftCutOffMean,
                                   aliveName = "alive",
                                   activityName = "activity",
                                   aAggregateCountName = "aAggregateCount",
                                   aDimerName =  "aDimer"
  )
  
  print(paste0("Logging"))
  sumAmyloidMonomers[i] <- sum(unlist(lapply(neurons, function(x) x$aMonomer))) / numberOfSynapses
  sumAmyloidDimers[i] <- sum(unlist(lapply(neurons, function(x) x$aDimer))) / numberOfSynapses
  sumAmyloidAggregatesCount[i] <- sum(unlist(lapply(neurons, function(x) x$aAggregateCount))) / numberOfSynapses
  sumAmyloidAggregatesSize[i] <- sum(unlist(lapply(neurons, function(x) x$aAggregateSum))) / numberOfSynapses 
  sumAmyoidPlaques[i] <- sum(unlist(lapply(neurons, function(x) x$aPlaque)))
  countAmyoidPlaques[i] <- sum(unlist(lapply(neurons, function(x) x$aPlaque > 0)))
  sumNFT[i] <- sum(unlist(lapply(neurons, function(x) x$nft))) / numberOfSynapses
  meanNFTSeedProbability[i] <- mean(unlist(lapply(neurons, function(x) x$nftSeedProbability)))
  meanNeuronActivity[i] <- mean(unlist(lapply(neurons, function(x) x$activity)))
  synapseLoss[i] <- mean(unlist(lapply(neurons, function(x) mean(x$activity == 0))))
  synapseActivity[i] <- mean(unlist(lapply(neurons, function(x) mean(x$activity))))
  alivePerc[i] <- mean(unlist(lapply(neurons, "[", "alive")))
  print(paste0("Percent Alive after current Step: ", round(alivePerc[i] * 100, 1), "%"))
  
  logModel(neurons = neurons, logPath =  logPath, step = i)
  logSampleNeuron(neurons = neurons, logPath =  logPath, t = i, index = c(1, 20, 40, 80, 100))
  
  #lapply(neurons, function(l){if(any(unlist(c(l$aMonomer, l$aDimer, l$aAggregateCount, l$aAggregateSum, l$aPlaque, l$nft, l$nftSeedProbability)) < 0)) stop(l$id)})
  #lapply(neurons, function(l){if(any(unlist(c(l$nftSeedProbability)) > 1)) stop(l$id)})
}
# lapply(neurons, function(l) {if(any(is.na(l$aMonomer))){print(paste("NA dected in aMonomer:")); print(l)}})
# lapply(neurons, function(l) {if(any(is.na(l$aDimer))){print(paste("NA dected in aDimer:")); print(l)}})
# lapply(neurons, function(l) {if(is.na(l$aDimerStack)){print(paste("NA dected in aDimerStack:")); print(l)}})
# lapply(neurons, function(l) {if(is.na(l$aMonomerStack)){print(paste("NA dected in aMonomerStack:")); print(l)}})

par(pch = 20)
#plot(sumAmyloidMonomers/(max(sumAmyloidMonomers)), col = "red", ylim = (c(0, max(sumAmyloidDimers, sumAmyloidMonomers, sumAmyloidAggregatesSize, na.rm = T))))
plot(sumAmyloidMonomers/(max(sumAmyloidMonomers)), col = "red", ylim = (c(0, 1)))
points(sumAmyloidDimers/(max(sumAmyloidDimers)), col = "blue")
points(sumAmyloidAggregatesSize/(max(sumAmyloidAggregatesSize)), col = "purple")
points(sumAmyloidAggregatesSize / sumAmyloidAggregatesCount / amyloidAggregateGeneration_aggregateMaxSize, col = "pink")
points((sumAmyoidPlaques / countAmyoidPlaques) / amyloidPlaqueGeneration_plaqueMaximumSize, col = "purple4")
points(alivePerc, col = "black")
points(synapseLoss, col = "blue4")
points(synapseActivity, col = "green")
points(sumNFT/max(sumNFT), col = "yellow")

plot(sumAmyloidAggregatesCount, col = "red", ylim = (c(0, max(sumAmyloidAggregatesCount, sumAmyloidAggregatesSize))))
points(sumAmyloidAggregatesSize, col = "blue")

plot(sumAmyloidAggregatesSize / sumAmyloidAggregatesCount / amyloidAggregateGeneration_aggregateMaxSize, col = "purple", ylim = c(0,1))
#points(1 - (sumAmyloidAggregatesCount / max(sumAmyloidAggregatesCount)), col = "purple4")
points((sumAmyoidPlaques / countAmyoidPlaques) / amyloidPlaqueGeneration_plaqueMaximumSize, col = "green")
points(alivePerc, col = "red")
points(synapseLoss, col = "blue")
points(synapseActivity, col = "yellow")

plot(sumAmyoidPlaques)
plot(countAmyoidPlaques / length(neurons), col = "red")
plot(sumNFT)
plot(meanNFTSeedProbability)
plot(sumNFT/synapseActivityUpdate_nftCutOffMean, ylim = c(0,1))
plot(alivePerc)
plot(sumAmyloidMonomers)
plot(sumAmyloidDimers)


# points(seq_along(sumAmyloidControle), sumAmyloidControle, col = "green")
# sumAmyloidControle <- sumAmyloid
# 
# sumAmyloidAD <- sumAmyloidControle
# 
# 
# neurons <- list.update(neurons, aMonomerStack = 0)
# resetAMonomers()
# 
# plot(sumAmyloidMonomers[1:500], col = "red", ylim = (c(0, max(sumAmyloidDimers[1:500], sumAmyloidMonomers[1:500]))))
# points(sumAmyloidDimers[1:500], col = "blue")