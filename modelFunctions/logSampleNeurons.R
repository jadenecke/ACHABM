logSampleNeuron <- function(neurons, t, index, logPath){
  if(t == 1){
    dir.create(file.path(logPath, "sampledNeurons"))
  }
  logPath <- file.path(logPath, "sampledNeurons")
  saveRDS(neurons[index], file = file.path(logPath, paste0("neurons_", t, ".rds")))
}
