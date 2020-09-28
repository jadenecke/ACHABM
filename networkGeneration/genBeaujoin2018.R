library(data.table)

if(file.exists("connectivityMatBeaujoin2018_shortNames.csv")){
  conMat <- as.matrix(fread("connectivityMatBeaujoin2018_shortNames.csv"), rownames = 1)
} else if (file.exists(file.path("..", "connectivityMatBeaujoin2018_shortNames.csv"))){
  conMat <- as.matrix(fread(file.path("..", "connectivityMatBeaujoin2018_shortNames.csv")), rownames = 1)
} else {
  stop("Beaujoin Mat not found, specify path")
}



conMat <- conMat[!grepl("FI", colnames(conMat), ignore.case = TRUE), !grepl("FI", colnames(conMat), ignore.case = TRUE)]
conMat <- conMat[!grepl("AL", colnames(conMat), ignore.case = TRUE), !grepl("AL", colnames(conMat), ignore.case = TRUE)]
conMat["H_EC",] <- conMat["H_EC",] + conMat["B_EC",]
conMat[,"H_EC"] <- conMat[,"H_EC"] + conMat[,"B_EC"]
conMat <- conMat[!grepl("B_EC", colnames(conMat), ignore.case = TRUE), !grepl("B_EC", colnames(conMat), ignore.case = TRUE)]

ratioHead <- conMat["H_CA1", "H_LAMO"] / conMat["H_CA23", "H_LAMO"]
ratioBody <- conMat["B_CA1", "B_LAMO"] / conMat["B_CA23", "B_LAMO"]
ratioTail <- conMat["T_CA1", "T_LAMO"] / conMat["T_CA23", "T_LAMO"]


# sum(conMat)
# sum(conMat[!grepl("LAMO", colnames(conMat), ignore.case = TRUE), !grepl("LAMO", colnames(conMat), ignore.case = TRUE)])
conMat["H_CA1",] <- conMat["H_CA1",] + round((ratioHead / (1+ratioHead)) * conMat["H_LAMO",])
conMat["H_CA23",] <- conMat["H_CA23",] + round((1 / (1+ratioHead)) * conMat["H_LAMO",])
conMat[,"H_CA1"] <- conMat[,"H_CA1"] + round((ratioHead / (1+ratioHead)) * conMat[,"H_LAMO"])
conMat[,"H_CA23"] <- conMat[,"H_CA23"] + round((1 / (1+ratioHead)) * conMat[,"H_LAMO"])

conMat["B_CA1",] <- conMat["B_CA1",] + round((ratioBody / (1+ratioBody)) * conMat["B_LAMO",])
conMat["B_CA23",] <- conMat["B_CA23",] + round((1 / (1+ratioBody)) * conMat["B_LAMO",])
conMat[,"B_CA1"] <- conMat[,"B_CA1"] + round((ratioBody / (1+ratioBody)) * conMat[,"B_LAMO"])
conMat[,"B_CA23"] <- conMat[,"B_CA23"] + round((1 / (1+ratioBody)) * conMat[,"B_LAMO"])

conMat["T_CA1",] <- conMat["T_CA1",] + round((ratioTail / (1+ratioTail)) * conMat["T_LAMO",])
conMat["T_CA23",] <- conMat["T_CA23",] + round((1 / (1+ratioTail)) * conMat["T_LAMO",])
conMat[,"T_CA1"] <- conMat[,"T_CA1"] + round((ratioTail / (1+ratioTail)) * conMat[,"T_LAMO"])
conMat[,"T_CA23"] <- conMat[,"T_CA23"] + round((1 / (1+ratioTail)) * conMat[,"T_LAMO"])

conMat <- conMat[!grepl("LAMO", colnames(conMat), ignore.case = TRUE), !grepl("LAMO", colnames(conMat), ignore.case = TRUE)]
diag(conMat) <- 0
conMatOriginal <- conMat
conMat <- round(sqrt(conMat))*25

# library(reshape2)
# library(ggplot2)
# conMatLong <- melt(conMat)
# ggplot(conMatLong) + 
#   geom_raster(aes(x = Var1, y = Var2, fill = value))



#EC size based on Price et al., 2001
#         H_EC      H_SUB     H_CA1   H_CA23   H_DG     B_SUB   B_CA1   B_CA23   B_DG     T_SUB   T_CA1   T_CA23  T_DG
size <- c(11850000, 1402170, 1223069, 313358, 2975444, 1123520, 1501601, 278920, 3240514, 144310, 2111441, 308864, 3244665)
sizeOriginal <- size
#size <- round(size / 17000)
size <- round(sqrt(size) / 10)
names(size) <- rownames(conMat)



# ###############
# #calculations (based on Cobb et al., 2013)
# #CA1: 
# round(4836111 * (11193/(11193 + 13742 + 19323)))#head
# round(4836111 * (13742/(11193 + 13742 + 19323)))#body
# round(4836111 * (19323/(11193 + 13742 + 19323)))#tail
# 
# #CA23:
# round(901142 * (16597/(16597 + 14773 + 16359)))#head
# round(901142 * (14773/(16597 + 14773 + 16359)))#body
# round(901142 * (16359/(16597 + 14773 + 16359)))#tail
# 
# #DG (+hilus)
# round((8028600 + 1432024) * ((11213 + 28928)/(((11213 + 28928) + (10778 + 32939) + (14405 + 29368)))))#head
# round((8028600 + 1432024) * ((10778 + 32939)/(((11213 + 28928) + (10778 + 32939) + (14405 + 29368)))))#body
# round((8028600 + 1432024) * ((14405 + 29368)/(((11213 + 28928) + (10778 + 32939) + (14405 + 29368)))))#tail
# 
# #SUB(based on Hrybousky et al., 2019 and West et al., 2000)
# round(2670000 * (133.6/(133.6 + 107.05 + 13.75)))#head
# round(2670000 * (107.05/(133.6 + 107.05 + 13.75)))#body
# round(2670000 * (13.75/(133.6 + 107.05 + 13.75)))#tail







