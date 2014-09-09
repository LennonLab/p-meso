################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment:  Microbial Community Indicator Species       #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/28                                                      #
#                                                                              #
################################################################################

pmeso.indval <- function(shared = " ", design = " ", level = "0.03"){

  # Setup Work Environment
  source("../functions/DiversityFunctions.r")
  require("vegan")
  require("labdsv")
  se <- function(x){sd(x)/sqrt(length(x))}

# Exp Design
  design <- read.delim(design)
  design$Treatment <- factor(design$Treatment,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Import Site by OTU Matrix
  pmeso_data <- t(read.otu(shared, level))

  # Remove Zero Sum OTUs
  pmeso_data <- pmeso_data[,!(colSums(abs(pmeso_data)) == 0)]

  # Calculate Presense Absence
  dataPA <- (pmeso_data > 0)*1

  # Calculating Relative Abundance
  dataREL <- pmeso_data
  for(i in 1:ncol(pmeso_data)){
    dataREL[,i] = pmeso_data[,i]/sum(pmeso_data[,i])
    }

  # Species Indicator Values
  pmeso.iva <- indval(t(dataREL), as.numeric(design$Treatment), numitr=1000)
  gr <- pmeso.iva$maxcls[pmeso.iva$pval <= 0.05]
  iv <- pmeso.iva$indcls[pmeso.iva$pval <= 0.05]
  pv <- pmeso.iva$pval[pmeso.iva$pval   <= 0.05]
  fr <- apply(t(dataREL)>0, 2, sum)[pmeso.iva$pval <= 0.05]
  fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
  fidg <- fidg[order(fidg$group, -fidg$indval), ]
  return(fidg)
  }
  
pmeso.indval.bio <- function(shared = " ", design = " ", 
  dummies = " ", level = "0.03"){

  # Setup Work Environment
  source("../functions/DiversityFunctions.r")
  require("vegan")
  require("labdsv")
  se <- function(x){sd(x)/sqrt(length(x))}

  # Exp Design
  design <- read.delim(design)
  design$Treatment <- factor(design$Treatment,
  levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Import Site by OTU Matrix
  pmeso_data <- t(read.otu(shared, level))
  
  dummy <- read.delim(dummies)
  dummy.trt1 <- dummy[dummy$cvr == "Rest",]
  dummy.trt1 <- dummy.trt1[dummy.trt1$mul == "Rest",]
    
  # Models of Interest
  samps1 <- as.character(dummy.trt1$Sample)

  # Import Site by OTU Matrix
  pmeso_data <- t(read.otu(shared, level))
  
  # Select only trt sites
  pmeso_data_samp <- subset(pmeso_data, select = samps1)

  # Remove Zero Sum OTUs
  pmeso_data <- pmeso_data[rowSums(abs(pmeso_data)) > 10,]
  pmeso_data_samp <- pmeso_data_samp[rowSums(abs(pmeso_data_samp)) > 10,]

  # Calculate Presense Absence
  dataPA <- (pmeso_data > 0)*1

  # Calculating Relative Abundance
  dataREL <- pmeso_data
  for(i in 1:ncol(pmeso_data)){
    dataREL[,i] = pmeso_data[,i]/sum(pmeso_data[,i])
    }
    
  dataREL2 <- pmeso_data_samp
  for(i in 1:ncol(pmeso_data_samp)){
    dataREL2[,i] = pmeso_data_samp[,i]/sum(pmeso_data_samp[,i])
    }

  # Species Indicator Values
  pmeso.iva <- indval(t(dataREL2), as.numeric(dummy.trt1$bio), numitr=1000)
  gr <- pmeso.iva$maxcls[pmeso.iva$pval <= 0.05]
  iv <- pmeso.iva$indcls[pmeso.iva$pval <= 0.05]
  pv <- pmeso.iva$pval[pmeso.iva$pval   <= 0.05]
  fr <- apply(t(dataREL2)>0, 2, sum)[pmeso.iva$pval <= 0.05]
  fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
  fidg <- fidg[order(fidg$group, -fidg$indval), ]
  return(fidg)
  }
