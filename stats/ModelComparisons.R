################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: Delta AIC/BIC Model Comparisons                 #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/05/05                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~./GitHub/p-meso/stats/")
require("car")
require("reshape")          

# Import Data 
# Eco
eco.data <- read.delim("../data/p-meso_eco.txt", header=T)
# Chem  
chem.data <- read.delim("../data/p-meso_chem.txt", header=T)
# Micro
bio.data <- read.delim("../data/p-meso_bio.txt", sep="\t", header=T)
# Exp Design
design <- read.delim("../data/design.txt")
  design$Treatment <- factor(design$Treatment, 
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

# Dummy Matrix
dummies <- read.delim("model_matrix.txt")

# Response Matrix
responses <- merge(merge(chem.data, bio.data, by="Sample"), 
  eco.data, by="Sample")

# Responses of Interest
response <- c("DOC.uM", "TN.uM", "TP.uM", "BP.uMC", "BR.uMO2", "BGE", "ChlA2", 
  "CladCope", "Resp_uM", "NPP_uM")
rownames(responses) <- responses$Sample
responses <- responses[as.character(dummies$Sample),]   # Removes Mix
responses <- responses[, response]                      # Only responses
  
# Models of Interest
#models <- c("y~trt+p+io+bio+com.b+com.g", 
# "y~trt+p+io+bio+com.b", "y~trt+p+io+bio+com.g", 
# "y~trt+p+io+bio", "y~trt+p+io+com.b", "y~trt+p+io+com.g", 
# "y~trt+p+io", "y~trt+p+bio", "y~trt+p+com.b", "y~trt+p+com.g",
# "y~p+io+bio+com.b+com.g", 
# "y~p+io+bio+com.b", "y~p+io+bio+com.g", 
# "y~p+io+bio", "y~p+io+com.b", "y~p+io+com.g", 
# "y~p+io", "y~p+bio", "y~p+com.b", "y~p+com.g",
# "y~trt+p", "y~trt+io", "y~trt+bio", "y~trt+com.g", "y~trt+com.b", "y~trt", "y~p", "y~io", "y~bio", "y~com.b", "y~com.g") 

models <- c("y~trt", "y~p", "y~io", "y~bio", "y~com.g")

# Output Matrix
out.aic <- matrix(data=NA, nrow=length(models),ncol=length(response))
row.names(out.aic) <- models
colnames(out.aic) <- response
out.bic <- matrix(data=NA, nrow=length(models),ncol=length(response))
row.names(out.bic) <- models
colnames(out.bic) <- response

# AIC & BIC Values
for (i in 1:length(response)){
  for (j in 1:length(models)){
    y = unlist(responses[response[i]])
    model.temp <- lm(as.formula(models[j]), data=dummies)
    out.aic[j,i] = AIC(model.temp)
    out.bic[j,i] = BIC(model.temp)
    }
  }
  
delta.aic <- matrix(data=NA, nrow=dim(out.aic)[1], ncol=dim(out.aic)[2])
row.names(delta.aic) <- row.names(out.aic)
colnames(delta.aic) <- colnames(out.aic)
for (i in 1:dim(out.aic)[2]){
  delta.aic[,i] <- round(out.aic[,i] - out.aic[1,i], 2)
  }
delta.aic

delta.aic < 2
delta.aic < 4


for (i in 1:dim(out.aic)[2]){
  delta.aic[,i] <- round(out.aic[,i] - out.aic[1,i], 2)
  }
delta.aic



