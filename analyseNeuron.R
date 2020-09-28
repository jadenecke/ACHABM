logFolder <- file.path("G:", "Uni", "SEDS MA Data", "logs")
logNumber <- "20200926_131402_brainAD"

logPath <- file.path(logFolder, logNumber, "sampledNeurons")
files <- list.files(logPath, full.names = T, pattern = ".rds")
library(rlist)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

#indexNumber <- 5

source(file.path(logFolder, logNumber, "modelParameter.R"))
non <- length(readRDS(files[1]))
dir.create(file.path(logPath, "neurons"))

for(i in seq(non)){
  neuron <- list()
  logStep <- as.numeric(lapply(strsplit(basename(files), "_|\\."), "[", 2))
  for(j in seq_along(logStep)){
    neuron[[logStep[j]]] <- readRDS(files[j])[[i]]
  }
  saveRDS(neuron, file.path(logPath, "neurons", paste0("neuron_", neuron[[1]]$areaName, "_", neuron[[1]]$id,".rds")))
}


#########################################
# 
# extractArray <- function(neuron, p){
#   a <- array(NA, dim = c( z = length(neuron), y = dim(neuron[[1]][[p]])[1], x = dim(neuron[[1]][[p]])[2]))
#   for(i in seq_along(neuron)){
#     a[i,,] <- neuron[[i]][[p]]
#   }
#   return(a)
# }
# 
# test <- extractArray(neuron, "aMonomer")
# testLong <- melt(test)
# testLong$value <- pmin(testLong$value, 20000)
# testLong$valScaled <- testLong$value/max(testLong$value)
# 
# testLongSubset <- testLong[testLong$Var1 %% (numberOfSteps/8) == 0,]
# testLongSubset$Var1 <- as.factor(testLongSubset$Var1)
# 
# 
# g1 <- ggplot(testLongSubset) + 
#   geom_raster(aes(y = Var3, x = Var2, fill = value))+
#   facet_wrap( ~ Var1, ncol = 8)+
#   labs(y = "Dendrites",
#        x = "Synapses") +
#   theme(panel.grid.minor = element_line(colour = "grey70", size = .5), legend.background = element_rect(fill="gray70", size=.5, linetype="dotted"))+
#   theme(strip.text = element_text(size = 11, face="bold"), strip.background = element_rect(color="gray80", fill="gray80", size=1.5, linetype="solid"))+
#   scale_y_continuous(expand = c(0, 0))+
#   scale_x_continuous(expand = c(0, 0))+
#   scale_fill_continuous(limits= c(0, max(testLong$value)), breaks = seq(0,max(testLong$value), by=max(testLong$value)/5))+
#   coord_fixed(ratio = max(testLongSubset$Var3)/max(testLongSubset$Var2), expand = TRUE, clip = "on")+
#   ggtitle("Before spread")+
#   theme(strip.text = element_text(size=8),
#         legend.key.height = unit(.4, "cm"),
#         # legend.box.background = element_rect(color="red", size=2),
#         legend.box.margin = margin(6, 1, 20, 1))
# 
# g2 <- ggplot(testLongSubset) + 
#   geom_raster(aes(y = Var3, x = Var2, fill = value))+
#   facet_wrap( ~ Var1, ncol = 8)+
#   labs(y = "Dendrites",
#        x = "Synapses") +
#   theme(panel.grid.minor = element_line(colour = "grey70", size = .5), legend.background = element_rect(fill="gray70", size=.5, linetype="dotted"))+
#   theme(strip.text = element_text(size = 11, face="bold"), strip.background = element_rect(color="gray80", fill="gray80", size=1.5, linetype="solid"))+
#   scale_y_continuous(expand = c(0, 0))+
#   scale_x_continuous(expand = c(0, 0))+
#   scale_fill_continuous(limits= c(0, max(testLong$value)), breaks = seq(0,max(testLong$value), by=max(testLong$value)/5))+
#   coord_fixed(ratio = max(testLongSubset$Var3)/max(testLongSubset$Var2), expand = TRUE, clip = "on")+
#   ggtitle("Before spread")+
#   theme(strip.text = element_text(size=8),
#         legend.key.height = unit(.4, "cm"),
#         # legend.box.background = element_rect(color="red", size=2),
#         legend.box.margin = margin(6, 1, 20, 1))
# 
# grid.arrange(g1, g2, ncol = 1)
