##########################################################
# 
##########################################################

amyloidMonomerCleaning <- function(neurons,
                                   mean,
                                   sd,
                                   upperBound,
                                   currentStep,
                                   declineFactor,
                                   declineFlag = TRUE,
                                   matName = "aMonomer",
                                   activityName = "activity" ){
  if(declineFlag){ #AD Model
    for(i in seq_along(neurons)){
      if(neurons[[i]][["alive"]]){
        neurons[[i]][[matName]] <- pmax(neurons[[i]][[matName]] - 
                                          pmin(round(
                                            rnorm(n = length(neurons[[i]][[matName]]), mean, sd) * 
                                              neurons[[i]][[activityName]] * 
                                              (1-rbeta(n = length(neurons[[i]][[matName]]), currentStep, declineFactor))
                                          ), upperBound), 0)
      }
    }
    return(neurons)
  } else { # Control Model
    for(i in seq_along(neurons)){
      if(neurons[[i]][["alive"]]){
        neurons[[i]][[matName]] <- pmax(neurons[[i]][[matName]] - 
                                          pmin(floor(
                                            rnorm(n = length(neurons[[i]][[matName]]), mean, sd) * 
                                              neurons[[i]][[activityName]]
                                          ), upperBound), 0)
      }
    }
    return(neurons)
  }
}


# amyloidMonomerCleaning_Cached <- function(neurons, # not faster
#                                    mean,
#                                    sd,
#                                    upperBound,
#                                    currentStep,
#                                    declineFactor,
#                                    declineFlag = TRUE,
#                                    matName = "aMonomer",
#                                    activityName = "activity" ){
#   
#   
#   if(
#     all(
#       unlist(
#         lapply(neurons,
#                function(l) is.null(l[[matName]])
#         )))){
#     stop(paste("Can't apply amyloidMonomerCleaning, either", matName, "or", activityName, "missing!"))
#   }
#   l <- max(unlist(lapply(neurons, function(l) length(l[["aMonomer"]])))) * 20
#   normPool <- rnorm(n = l, mean, sd)
#   betaPool <- 1-rbeta(n = l, currentStep, declineFactor)
#   
#   if(declineFlag){ #AD Model
#     for(i in seq_along(neurons)){
#       neurons[[i]][[matName]] <- pmax(neurons[[i]][[matName]] - 
#                                         pmin(floor(
#                                           sample(normPool, length(neurons[[i]][[matName]])) * 
#                                             neurons[[i]][[activityName]] * 
#                                             sample(betaPool, length(neurons[[i]][[matName]]))
#                                         ), upperBound), 0)
#     }
#     return(neurons)
#   } else { # Control Model
#     for(i in seq_along(neurons)){
#       neurons[[i]][[matName]] <- pmax(neurons[[i]][[matName]] - 
#                                         pmin(floor(
#                                           sample(normPool, length(neurons[[i]][[matName]])) * 
#                                             neurons[[i]][[activityName]]
#                                         ), upperBound), 0)
#     }
#     return(neurons)
#   }
#   
# }
