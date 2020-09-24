##########################################################
# This script defines all the hard-coded hyper-parameter
##########################################################

#### Model ####
numberOfSteps <- 2500 #[1 - Inf)

#### Network Generation ####

# Parameter (very - small)
numberOfCluster <- 2 #number of areas
clusterMinSize <- 55
clusterMaxSize <- 70
synapsesPerDendriteMin <- 8
synapsesPerDendriteMax <- 11
dendritesPerNeuronMin <- 4
dendritesPerNeuronMax <- 5

# # Parameter (small)
# numberOfCluster <- 3 #number of areas
# clusterMinSize <- 90
# clusterMaxSize <- 120
# synapsesPerDendriteMin <- 8
# synapsesPerDendriteMax <- 15
# dendritesPerNeuronMin <- 4
# dendritesPerNeuronMax <- 6

# # Parameter (medium-small)
# numberOfCluster <- 3 #number of areas
# clusterMinSize <- 350
# clusterMaxSize <- 500
# synapsesPerDendriteMin <- 20
# synapsesPerDendriteMax <- 40
# dendritesPerNeuronMin <- 3
# dendritesPerNeuronMax <- 7

# # Parameter (medium-large)
# n <- 3 #number of areas
# clusterMinSize <- 1000
# clusterMaxSize <- 1500
# synapsesPerDendriteMin <- 30
# synapsesPerDendriteMax <- 100
# dendritesPerNeuronMin <- 5
# dendritesPerNeuronMax <- 10
# # 
# # # Parameter (large)
# n <- 7 #number of areas
# clusterMinSize <- 3000
# clusterMaxSize <- 4000
# synapsesPerDendriteMin <- 250
# synapsesPerDendriteMax <- 500
# dendritesPerNeuronMin <- 5
# dendritesPerNeuronMax <- 9

#### Amyloid Monomer Generation ####
amyloidMonomerGeneration_mu <- 3000 # [1 - Inf) mean number of amyloid monomer generation per time step per synapse
amyloidMonomerGeneration_delta <- 1000# [1 - Inf) standard deviation of amyloid monomer generation, 66% of synapses generate mu+-delta monomers per time step
amyloidMonomerGeneration_upperBound <- amyloidMonomerGeneration_mu + (5 * amyloidMonomerGeneration_delta) # [1 - Inf) upper boundary for new amyloid monomer addition per time step

#### Amyloid Monomer Clearance  ####
amyloidMonomerCleaning_mu <- amyloidMonomerGeneration_mu * 1.12 # [1 - Inf) mean number of amyloid monomer generation per time step per synapse
amyloidMonomerCleaning_delta <- amyloidMonomerGeneration_delta # [1 - Inf) standard deviation of amyloid monomer generation, 66% of synapses generate mu+-delta monomers per time step
amyloidMonomerCleaning_upperBound <- amyloidMonomerGeneration_upperBound # [1 - Inf)
amyloidMonomerCleaning_declineFactor <- numberOfSteps * 1.5 # (0 - 1) results in > 50% decline, (1 - Inf) results in < 50% decline,  results in ~50% reduction at the last time step, if numberOfSteps * 2 = ~ 33%, numberOfSteps / 2 = ~ 66%
# plot(1:numberOfSteps, sapply(1:numberOfSteps, function(x) 1-rbeta(1, x, numberOfSteps * 2)), ylim = c(0,1)); abline(.6666, 0)
# plot(1-as.vector(sapply(1:100, function(i) rbeta(100, i, 100*.5))))

#### Amyloid Monomer Diffusion - Intra Dendrite ####
amyloidMonomerDiffusionIntraDendrite_amountOfDiffusion <- .1 #[0-1) percent of monomers that will defuse away from the neuron (amount is before randomization, s.t the actual amount might vary)
amyloidMonomerDiffusionIntraDendrite_range <- 4 #[1 - Inf) range of spread to the left and right of the origin
amyloidMonomerDiffusionIntraDendrite_spreadSD <- .6 #[1 - Inf) #SD of normal distribution of the spread. higher values mean more even spread
amyloidMonomerDiffusionIntraDendrite_spreadMaxMultiplyer <- 4 #[1 - Inf) # factor introducing randomness: multiplying factor for each synapse spread with the range of 1/value - value, higher values mean more variation in the spread variation

#### Amyloid Monomer Diffusion - Inter Dendrite ####
amyloidMonomerDiffusionInterDendrite_amountOfDiffusion <- .02 #[0-1) percent of monomers that will defuse away from the neuron (amount is before randomization, s.t the actual amount might vary)
amyloidMonomerDiffusionInterDendrite_range <- 2 #[1 - Inf) range of spread to the left and right of the origin
amyloidMonomerDiffusionInterDendrite_spreadSD <- .6 #[1 - Inf) #SD of normal distribution of the spread. higher values mean more even spread
amyloidMonomerDiffusionInterDendrite_spreadMaxMultiplyer <- 4 #[1 - Inf) # factor introducing randomness: multiplying factor for each synapse spread with the range of 1/value - value, higher values mean more variation in the spread variation

#### Amyloid Monomer Diffusion - Inter Neuron ####
# The monomers are drawn from every snyapse based on a binomial distribution.
# The probability of the binomial distribution is calculated based on a logistic curve, s.t. higher amyloid loads lead to a higher probability of amyloid spread.
amyloidMonomerDiffusionInterNeuron_maximumSpreadProbability <- .05# [0:1) precent of monomers that can be spread if load is maximal
amyloidMonomerDiffusionInterNeuron_spreadDependencyCurveSteepness <- .001 # [0:1] determines the steepness of the logistic curve: values of ~1 or larger result in a sharp cutoff arround the inflection point with lower amyloid load resulting in 0 spread and higher amyloid load resulting in the maximum spread. Values around 0.01 - 0.001 result in a smooth increase within an amyloid range of 0 : 10,000. Also highly depends on the Inflection Point, especially for low amyloid loads, so check the curve carefully
amyloidMonomerDiffusionInterNeuron_spreadDependencyCurveInflectionPoint <- 6000# [1:Inf)
# curve(.2/(1 + exp(-.001*(x - 2000))), 1, 4000)
# curve(L/(1 + exp(-k*(x - x_zero))), 1, 10000) # L = Maxima; k = steepness; x_zero = inflection Point (Wendepunkt)

#### Amyloid Dimer Generation ####
amyloidDimerGeneration_maximumPercentTransform <- .1# [0:1) precent of monomers that can be converted if load is maximal
amyloidDimerGeneration_probCurveSteepness <- .001 # [0:1] determines the steepness of the logistic curve: values of ~1 or larger result in a sharp cutoff arround the inflection point with lower amyloid load resulting in 0 spread and higher amyloid load resulting in the maximum spread. Values around 0.01 - 0.001 result in a smooth increase within an amyloid range of 0 : 10,000. Also highly depends on the Inflection Point, especially for low amyloid loads, so check the curve carefully
amyloidDimerGeneration_probCurveInflectionPoint <- 10000 # [1:Inf)
#curve(.1/(1 + exp(-.001*(x - 3000))), 1, 6000)

######################################################

#### Amyloid Dimer Diffusion - Intra Dendrite ####
amyloidDimerDiffusionIntraDendrite_amountOfDiffusion <- .02 #[0-1) percent of Dimer that will defuse away from the neuron (amount is before randomization, s.t the actual amount might vary)
amyloidDimerDiffusionIntraDendrite_range <- 2 #[1 - Inf) range of spread to the left and right of the origin
amyloidDimerDiffusionIntraDendrite_spreadSD <- .6 #[1 - Inf) #SD of normal distribution of the spread. higher values mean more even spread
amyloidDimerDiffusionIntraDendrite_spreadMaxMultiplyer <- 2 #[1 - Inf) # factor introducing randomness: multiplying factor for each synapse spread with the range of 1/value - value, higher values mean more variation in the spread variation

#### Amyloid Dimer Diffusion - Inter Dendrite ####
amyloidDimerDiffusionInterDendrite_amountOfDiffusion <- .02 #[0-1) percent of Dimer that will defuse away from the neuron (amount is before randomization, s.t the actual amount might vary)
amyloidDimerDiffusionInterDendrite_range <- 1 #[1 - Inf) range of spread to the left and right of the origin
amyloidDimerDiffusionInterDendrite_spreadSD <- .6 #[1 - Inf) #SD of normal distribution of the spread. higher values mean more even spread
amyloidDimerDiffusionInterDendrite_spreadMaxMultiplyer <- 2 #[1 - Inf) # factor introducing randomness: multiplying factor for each synapse spread with the range of 1/value - value, higher values mean more variation in the spread variation

#### Amyloid Dimer Diffusion - Inter Neuron ####
# The monomers are drawn from every snyapse based on a binomial distribution.
# The probability of the binomial distribution is calculated based on a logistic curve, s.t. higher amyloid loads lead to a higher probability of amyloid spread.
amyloidDimerDiffusionInterNeuron_maximumSpreadProbability <- .1# [0:1) precent of Dimer that can be spread if load is maximal
amyloidDimerDiffusionInterNeuron_spreadDependencyCurveSteepness <- .001 # [0:1] determines the steepness of the logistic curve: values of ~1 or larger result in a sharp cutoff arround the inflection point with lower amyloid load resulting in 0 spread and higher amyloid load resulting in the maximum spread. Values around 0.01 - 0.001 result in a smooth increase within an amyloid range of 0 : 10,000. Also highly depends on the Inflection Point, especially for low amyloid loads, so check the curve carefully
amyloidDimerDiffusionInterNeuron_spreadDependencyCurveInflectionPoint <- 1000# [1:Inf)
# curve(.2/(1 + exp(-.001*(x - 2000))), 1, 4000)
# curve(L/(1 + exp(-k*(x - x_zero))), 1, 10000) # L = Maxima; k = steepness; x_zero =     inflection Point (Wendepunkt)

#### Amyloid Dimer Dissaggregation ####
amyloidDimerDisaggregation_DisaggregationProbability <- .10 # flat probability of a binomial distribution

#### Amyloid Aggregate Generation ####
amyloidAggregateGeneration_aggregateMaxSize <- 24
amyloidAggregateGeneration_seedingProbabiltyMax <- .03
amyloidAggregateGeneration_seedingProbabilityCurveSteepness <- .0002
amyloidAggregateGeneration_seedingProbabilityCurveInflectionPoint <- 12000
#curve(.004/(1 + exp(-.002*(x - 3000))), 1, 10000)
amyloidAggregateGeneration_aggregateGrowthProbabiltyMax <- .2
amyloidAggregateGeneration_aggregateGrowthProbabilityCurveSteepness <- .0005
amyloidAggregateGeneration_aggregateGrowthDelay <- 2000
amyloidAggregateGeneration_maxDecline <- .1 # [0 - 1], real max Decline seems to be maxDecline/2 as the inflections point for the highest declining curve is at b
amyloidAggregateGeneration_stability <- 4 # (0 - Inf) stability, hard to describe, but the range of 2-6 seems reasonable, at higher or lower values, the dependencie on the percent saturation becomes neglectable 

#### Amyloid Plaque Generation ####
amyloidPlaqueGeneration_plaqueMaximumSize <- 10000000 # [1 , Inf] maximum number of monomer units a plaque can hold
amyloidPlaqueGeneration_plaqueSeedProbabilityCurveMax <- .0005 # (0 - 1], real max Probability seems to be at x/2
amyloidPlaqueGeneration_plaqueSeedProbabilityCurveSteepness <- .01  #(0 - Inf) seems to relate to the number of aggregates at which the probability as roughly passed 50% of the max/2 probability
amyloidPlaqueGeneration_plaqueSeedProbabilityInflectionPoint <- 700 # inflecton point: dimers + aggregate count
amyloidPlaqueGeneration_plaquePullIntraDendrite <- 11 #[0 , Inf] range of intra dendrite pull from other synapses into the plaque
amyloidPlaqueGeneration_plaquePullInterDendrite <- 5 #[0 , Inf] range of inter dendrite pull from other synapses into the plaque
amyloidPlaqueGeneration_plaquePullMaxProb <- 0.5 #[0 , 1] maximum percent of monomers/dimer/aggregate which can be pulled in one step
amyloidPlaqueGeneration_AggregatesizePullRelation <- .5 #(0 , 1) values close to 0 indicate only large Aggregates are used to fill Plaques
amyloidPlaqueGeneration_softLimit <- 200 # [0, Inf] soft limit (maximumSize - softLimit) at which a plaque already stops growing in favor of computational performance

#### NFT Generation ####
nftGeneration_nftSeedProbability <- .0008 # Initial value for the NFT seed probability. The size from which the generated NFT is drawn is dependent on the aAggregateCount and nft itself
nftGeneration_nftAcceleration <- .02 #[0, Inf] formula is: sqrt(nft + aggregateSum) / Acceleration. Values < 1 are a speedup, values > 1 is slowdown
nftGeneration_maxNftGrowth <- 2 # maximum percent of nft growth per step of the nftCutOffMean (but has to be multiplies with the seed proability to yield actual value)
nftGeneration_nftFlatClearance <- 2 #[0, Inf] flat number of NFT that can be maximally cleaned in one step, with nftFlatClearance/2 * activity being the expected value (beta distribution)

#### NFT Seed Probability Inter Neuron Spread ####
nftSeedProbabilitySpread_lag <- 300000000 # [0, Inf] Count of NFT at which the probability rise is .5
nftSeedProbabilitySpread_maxSeedProb <- .003 #[0, 1] maximum nft seed probability which can be reached by a neuron
nftSeedProbabilitySpread_seedProbMaxRise <- .00005 # [0, 1] maximum nftSeedProbability rise in one step 
#plot(rbeta(10000, 10000, nftSeedProbabilitySpread_lag))
#plot(unlist(lapply(1:synapseActivityUpdate_nftCutOffMean, function(n){75 * rbeta(1, n , nftSeedProbabilitySpread_lag)}))) # probability addition per step depending on number of NFT


#### Synapse / Neurons Activity Update #### 
synapseActivityUpdate_declineStability <- 4000 #shape1 parameter of the beta distribution 
synapseActivityUpdate_declineStart <- 20000000  #shape2 parameter of the beta distribution 
synapseActivityUpdate_activityCutOff <- .1 #cut-off value at which the neuron activity is switched to zero and the neuron is also considered dead
synapseActivityUpdate_nftCutOffMean <- 20000 #value at which a neuron is considered dead (is calculated as the mean value of the nft matrix). At 20.000 there can be no more then 107374 synapses per neuron, otherwise the model does not work (Integer limit: 2^31 / 20000)
#plot(unlist(lapply(1:20000, function(x) rbeta(1, declineStability/sqrt(x/10), x^2 / declineStart))), ylim = c(0,1));abline(h = .2)


####  ####
