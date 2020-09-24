##########################################################
# 
##########################################################

logParamter <- function(logPath){
  #copy Model parameter
  file.copy("modelParameter.R", file.path(logPath, "modelParameter.R"))
  return(NULL)
}
