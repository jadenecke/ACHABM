#logFolder <- file.path("results")
logFolder <- file.path("G:", "Uni", "SEDS MA Data", "logs")
logNumber <- "20200927_024945_brainADIntervention1"

logPath <- file.path(logFolder, logNumber)
files <- list.files(logPath, full.names = T, pattern = ".rds")
library(rlist)

source(file.path(logPath, "modelParameter.R"))
#assign("x", 5)
logs <- list()
logStep <- as.numeric(lapply(strsplit(basename(files), "_|\\."), "[", 2))
for(i in seq_along(logStep)){
  logs[[logStep[i]]] <- readRDS(files[i])
}

getIndicatorPerNeuron <- function(logs, ind, func){
  areas <- unique(unlist(lapply(logs[[1]][["neurons"]], "[", "areaName")))
  df <- matrix(NA, nrow = length(logs), ncol = length(areas))
  colnames(df) <- areas
  pb = txtProgressBar(min = 0, max = length(logs), initial = 0, style = 3) 
  for(i in seq_along(logs)){
    setTxtProgressBar(pb,i)
    for(j in seq_along(areas)){
      df[i,j] <- func(unlist(lapply(logs[[i]][["neurons"]], function(n) if(n$areaName == areas[j]){n[[ind]]} else {return(NULL)})))
    }
  }
  return(df)
}

getIndicatorPerNeuronSizeStandardized <- function(logs, ind, func){
  areas <- unique(unlist(lapply(logs[[1]][["neurons"]], "[", "areaName")))
  df <- matrix(NA, nrow = length(logs), ncol = length(areas))
  colnames(df) <- areas
  pb = txtProgressBar(min = 0, max = length(logs), initial = 0, style = 3) 
  for(i in seq_along(logs)){
    setTxtProgressBar(pb,i)
    for(j in seq_along(areas)){
      df[i,j] <- func(unlist(lapply(logs[[i]][["neurons"]], function(n) if(n$areaName == areas[j]){n[[ind]]} else {return(NULL)})))
    }
  }
  return(df)
}

getIndicatorPerSynapse <- function(logs, ind){
  areas <- unique(unlist(lapply(logs[[1]][["neurons"]], "[", "areaName")))
  df <- matrix(NA, nrow = length(logs), ncol = length(areas))
  colnames(df) <- areas
  pb = txtProgressBar(min = 0, max = length(logs), initial = 0, style = 3) 
  for(i in seq_along(logs)){
    setTxtProgressBar(pb,i)
    for(j in seq_along(areas)){
      df[i,j] <- mean(unlist(lapply(logs[[i]][["neurons"]], function(n) if(n$areaName == areas[j]){n[[ind]] / n[["numberOfsynapses"]]} else {return(NULL)})))
    }
  }
  return(df)
}



results <- list(
  perNeuron = list(
    numberOfsynapses = getIndicatorPerNeuron(logs, "numberOfsynapses", mean),
    aMonomerMean = getIndicatorPerNeuron(logs, "aMonomerSum", mean),
    aDimerMean = getIndicatorPerNeuron(logs, "aDimerSum", mean),
    aAggregateCountMean = getIndicatorPerNeuron(logs, "aAggregateCountSum", mean),
    aAggregateSumMean = getIndicatorPerNeuron(logs, "aAggregateSumSum", mean),
    aPlaqueNumberMean = getIndicatorPerNeuron(logs, "aPlaqueNumber", mean),
    aPlaqueSizeMean = getIndicatorPerNeuron(logs, "aPlaqueSum", mean),
    nftSeedProbabilityMean = getIndicatorPerNeuron(logs, "nftSeedProbability", mean),
    nftSumMean = getIndicatorPerNeuron(logs, "nftSum", mean),
    percentAlive = getIndicatorPerNeuron(logs, "alive", mean),
    meanActivity = getIndicatorPerNeuron(logs, "activity", mean),
    synapseAlivePercentMean = getIndicatorPerNeuron(logs, "synapseAlivePercent", mean)
  ),
  perSynapse = list(
    aMonomerMean = getIndicatorPerSynapse(logs, "aMonomerSum"),
    aDimerMean = getIndicatorPerSynapse(logs, "aDimerSum"),
    aAggregateCountMean = getIndicatorPerSynapse(logs, "aAggregateCountSum"),
    aAggregateSumMean = getIndicatorPerSynapse(logs, "aAggregateSumSum"),
    nftSumMean = getIndicatorPerSynapse(logs, "nftSum")
  )
)

fileName <- paste0(logNumber, ".rds")
if(!file.exists(file.path("results", fileName))){saveRDS(results, file.path("results", fileName))}else{stop("File already exists!")}
