################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment:  Microbial Community PERMANOVA               #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/27                                                      #
#                                                                              #
################################################################################

pmeso.permanova <- function(shared = " ", design = " ", level = "0.03"){

  # Setup Work Environment
  source("../functions/DiversityFunctions.r")
  require("vegan")
  se <- function(x){sd(x)/sqrt(length(x))}
  
  # Exp Design
  design <- read.delim(design)
  design$Treatment <- factor(design$Treatment,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Import Site by OTU Matrix
  pmeso_data <- t(read.otu(shared, level))

  # Remove Zero Sum OTUs
  pmeso_data <- pmeso_data[,!(colSums(abs(pmeso_data)) ==0)]

  # Calculate Presense Absence
  dataPA <- (pmeso_data > 0)*1

  # Calculating Relative Abundance
  dataREL <- pmeso_data
  for(i in 1:ncol(pmeso_data)){
    dataREL[,i] = pmeso_data[,i]/sum(pmeso_data[,i])
    }

  # Create Distance Matrix
  samplePA.dist <- vegdist(decostand(t(dataPA),method="log"),method="bray")
  sampleREL.dist <- vegdist(decostand(t(dataREL),method="log"),method="bray")
  
  # Adonis (PERMANOVA)
  # Adonis runs a PERMANOVA (Created by Marti J. Anderson)
  # this is very similar to ANOVA but for multivariate data
  # You can make very complex experimental designs with it
  # The default distance measure is bray-curtis, but other measures
  # (Chao, Jaccard, Euclidean) can be used when specified
  # I am using a loop to run Adonis for each design model
  
  Adonis <- adonis(sampleREL.dist ~ design$Treatment, method="bray",
    permutations=9999)
    return(Adonis)
  }
  
pmeso.multi.permanova <- function(shared = " ", design = " ", 
  dummies = " ", level = "0.03"){

  # Temp
  # shared = "../mothur/output/tanks.bac.final.shared"
  # design = "../data/design.txt"
  # dummies = "../stats/p-meso_dummies.txt"
  # level = "0.03"

  # Setup Work Environment
  source("../functions/DiversityFunctions.r")
  require("vegan")
  se <- function(x){sd(x)/sqrt(length(x))}
  
  # Exp Design
  design <- read.delim(design)
  design$Treatment <- factor(design$Treatment,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))
    
  dummy <- read.delim(dummies)
  dummy.trt1 <- dummy[dummy$cvr == "Rest",]
  dummy.trt2 <- droplevels(dummy.trt1[dummy.trt1$mul == "Rest",])
    
  # Models of Interest
  models <- colnames(dummy.trt2)[2:11]
  samps1 <- as.character(dummy.trt1$Sample)
  samps2 <- as.character(dummy.trt2$Sample)

  # Import Site by OTU Matrix
  pmeso_data <- t(read.otu(shared, level))
  
  # Select only trt sites
  pmeso_data_samp <- subset(pmeso_data, select = samps1)
  pmeso_data_trt <- subset(pmeso_data, select = samps2)

  # Remove Zero Sum OTUs
  pmeso_data_samp <- pmeso_data_samp[,!(colSums(abs(pmeso_data_samp)) == 0)]
  pmeso_data_trt <- pmeso_data_trt[,!(colSums(abs(pmeso_data_trt)) == 0)]

  # Calculate Presense Absence
  dataPA <- (pmeso_data_trt > 0)*1

  # Calculating Relative Abundance
  dataREL1 <- pmeso_data_samp
  for(i in 1:ncol(pmeso_data_samp)){
    dataREL1[,i] = pmeso_data_samp[,i]/sum(pmeso_data_samp[,i])
    }
    
  dataREL2 <- pmeso_data_trt
  for(i in 1:ncol(pmeso_data_trt)){
    dataREL2[,i] = pmeso_data_trt[,i]/sum(pmeso_data_trt[,i])
    }
    
  dataREL3 <- pmeso_data
  for(i in 1:ncol(pmeso_data)){
    dataREL3[,i] = pmeso_data[,i]/sum(pmeso_data[,i])
    }
                                             
  # Create Distance Matrix
  samplePA.dist <- vegdist(decostand(t(dataPA),method="log"),method="bray")
  sampleREL1.dist <- vegdist(decostand(t(dataREL1),method="log"),method="bray")
  sampleREL2.dist <- vegdist(decostand(t(dataREL2),method="log"),method="bray")
  sampleREL3.dist <- vegdist(decostand(t(dataREL3),method="log"),method="bray")
  
  # Adonis (PERMANOVA)
  # Adonis runs a PERMANOVA (Created by Marti J. Anderson)
  # this is very similar to ANOVA but for multivariate data
  # You can make very complex experimental designs with it
  # The default distance measure is bray-curtis, but other measures
  # (Chao, Jaccard, Euclidean) can be used when specified
  # I am using a loop to run Adonis for each design model
  
  adonis.t <- adonis(sampleREL2.dist ~ dummy.trt2$Treatment, method="bray", permutations=999)
  adonis.p <- adonis(sampleREL3.dist ~ dummy$cvr, method="bray", permutations=999)
  adonis.i <- adonis(sampleREL2.dist ~ dummy.trt2$io, method="bray", permutations=999)
  adonis.b <- adonis(sampleREL2.dist ~ dummy.trt2$bio, method="bray", permutations=999)
  adonis.c <- adonis(sampleREL2.dist ~ dummy.trt2$com.g, method="bray", permutations=999)
  adonis.m <- adonis(sampleREL1.dist ~ dummy.trt1$mul, method="bray", permutations=999)
  
  return(list(adonis.t, adonis.p, adonis.i, adonis.b, adonis.c, adonis.m)) 
  }