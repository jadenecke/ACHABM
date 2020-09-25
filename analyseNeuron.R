logFolder <- file.path("G:", "Uni", "SEDS MA Data", "logs")
logNumber <- "20200924_162544"

logPath <- file.path(logFolder, logNumber, "sampledNeurons")
files <- list.files(logPath, full.names = T, pattern = ".rds")
library(rlist)

indexNumber <- 1

source(file.path(logFolder, logNumber, "modelParameter.R"))
#assign("x", 5)
neuron <- list()
logStep <- as.numeric(lapply(strsplit(basename(files), "_|\\."), "[", 2))
for(i in seq_along(logStep)){
  neuron[[logStep[i]]] <- readRDS(files[i])[[indexNumber]]
}


extractArray <- function(neuron, p){
  a <- array(NA, dim = c( z = length(neuron), y = dim(neuron[[1]][[p]])[1], x = dim(neuron[[1]][[p]])[2]))
  for(i in seq_along(neuron)){
    a[i,,] <- neuron[[i]][[p]]
  }
  return(a)
}

test <- extractArray(neuron, "activity")
testLong <- melt(test)
testLong$value <- pmin(testLong$value, 20000)
testLong$valScaled <- testLong$value/max(testLong$value)

testLongSubset <- testLong[testLong$Var1 %% 400 == 0,]
testLongSubset$Var1 <- as.factor(testLongSubset$Var1)


ggplot(testLongSubset) + 
  geom_raster(aes(y = Var3, x = Var2, fill = value))+
  facet_wrap( ~ Var1, ncol = 5)+
  labs(y = "Dendrites",
       x = "Synapses") +
  theme(panel.grid.minor = element_line(colour = "grey70", size = .5), legend.background = element_rect(fill="gray70", size=.5, linetype="dotted"))+
  theme(strip.text = element_text(size = 11, face="bold"), strip.background = element_rect(color="gray80", fill="gray80", size=1.5, linetype="solid"))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  scale_fill_continuous(limits= c(0, max(testLong$value)), breaks = seq(0,max(testLong$value), by=max(testLong$value)/5))+
  coord_fixed(ratio = .45, expand = TRUE, clip = "on")+
  ggtitle("Before spread")
