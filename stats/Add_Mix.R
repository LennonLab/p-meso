################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: Additive Response, deviance models              #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/04/29                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~./GitHub/p-meso/stats/")

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

# Response Matrix
responses <- merge(merge(chem.data, bio.data, by="Sample"), 
  eco.data, by="Sample")
    
# Responses of Interest
response <- c("DOC.uM", "TN.uM", "TP.uM", "BP.uMC", "BR.uMO2", "BGE", "ChlA2", 
  "CladCope", "Resp_uM", "NPP_uM")
  
responses <- responses[, response]

# Calculate Deviance (Loreau 1998)

dev.tests <- matrix(NA,length(response),3)
colnames(dev.tests) <- c("Response", "Avg Deviance", "P-value")
dev.tests[,1] <- response
for (i in 1:length(response)){
  y.e <- (0.25*(tapply(responses[,i], design$Treatment, mean, simplify=F)$SRP) + 
          0.25*(tapply(responses[,i], design$Treatment, mean, simplify=F)$AEP) +
          0.25*(tapply(responses[,i], design$Treatment, mean, simplify=F)$ATP) +
          0.25*(tapply(responses[,i], design$Treatment, mean, simplify=F)$Phyt))
  y.o <- responses[,i][design$Treatment == "Mix"]
  dev <- rep(NA, length(y.o))
  for (j in 1:length(y.o)){
    dev[j] = (y.o[j] - y.e)/(y.e)
    }
  dev.tests[i,2] <- mean(dev, na.rm=T)
  dev.tests[i,3] <- t.test(dev)$p.value
  }
  
  
################################################################################
# Diversity

source("../functions/p-meso.diversity.R")

# Heterotrophs
heteros <- pmeso.diversity.calc(shared = "../mothur/output/tanks.bac.final.shared",
  design = "../data/design.txt", level = "0.03", size = 10000)

# Phototrophs
cyanos <- pmeso.diversity.calc(shared = "../mothur/output/tanks.photo.final.shared",
  design = "../data/design.txt", level = "0.03", size = 10000)

# Eukaryotes
euks <- pmeso.diversity.calc(shared = "../mothur/output/tanks.euk.final.shared",
  design = "../data/design.txt", level = "unique", size = 100)

# Calculate Deviance (Loreau 1998)

div.dev.tests <- matrix(NA,3,3)
colnames(div.dev.tests) <- c("Response", "Avg Deviance", "P-value")
div.dev.tests[,1] <- c("heteros", "cyanos", "euks")


  y.e <- (0.25*(tapply(heteros$Richness$mean, design$Treatment, mean, simplify=F)$SRP) + 
          0.25*(tapply(heteros$Richness$mean, design$Treatment, mean, simplify=F)$AEP) +
          0.25*(tapply(heteros$Richness$mean, design$Treatment, mean, simplify=F)$ATP) +
          0.25*(tapply(heteros$Richness$mean, design$Treatment, mean, simplify=F)$Phyt))
  y.o <- heteros$Richness$mean[design$Treatment == "Mix"]
  dev <- rep(NA, length(y.o))
  for (j in 1:length(y.o)){
    dev[j] = (y.o[j] - y.e)/(y.e)
    }
  div.dev.tests[1,2] <- mean(dev, na.rm=T)
  div.dev.tests[1,3] <- t.test(dev)$p.value

  
  y.e <- (0.25*(tapply(cyanos$Richness$mean, design$Treatment, mean, simplify=F)$SRP) + 
          0.25*(tapply(cyanos$Richness$mean, design$Treatment, mean, simplify=F)$AEP) +
          0.25*(tapply(cyanos$Richness$mean, design$Treatment, mean, simplify=F)$ATP) +
          0.25*(tapply(cyanos$Richness$mean, design$Treatment, mean, simplify=F)$Phyt))
  y.o <- cyanos$Richness$mean[design$Treatment == "Mix"]
  dev <- rep(NA, length(y.o))
  for (j in 1:length(y.o)){
    dev[j] = (y.o[j] - y.e)/(y.e)
    }
  div.dev.tests[2,2] <- mean(dev, na.rm=T)
  div.dev.tests[2,3] <- t.test(dev)$p.value
  

  y.e <- (0.25*(tapply(euks$Richness$mean, design$Treatment, mean, simplify=F)$SRP) + 
          0.25*(tapply(euks$Richness$mean, design$Treatment, mean, simplify=F)$AEP) +
          0.25*(tapply(euks$Richness$mean, design$Treatment, mean, simplify=F)$ATP) +
          0.25*(tapply(euks$Richness$mean, design$Treatment, mean, simplify=F)$Phyt))
  y.o <- euks$Richness$mean[design$Treatment == "Mix"]
  dev <- rep(NA, length(y.o))
  for (j in 1:length(y.o)){
    dev[j] = (y.o[j] - y.e)/(y.e)
    }
  div.dev.tests[3,2] <- mean(dev, na.rm=T)
  div.dev.tests[3,3] <- t.test(dev)$p.value
  
  
  
  