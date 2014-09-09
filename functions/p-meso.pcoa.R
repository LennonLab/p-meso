################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment:  Microbial Community PcoA and PERMANOVA      #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/04/21                                                      #
#                                                                              #
################################################################################

pmeso.pcoa <- function(shared = " ", design = " ", plot.title = "test", 
  level = "0.03"){

  # Setup Work Environment
  source("../functions/DiversityFunctions.r")  
  require("vegan")
  se <- function(x){sd(x)/sqrt(length(x))}                                   

  # Import Site by OTU Matrix
  pmeso_data <- t(read.otu(shared, level))             
  design <- read.delim(design, header=T, row.names=1)
  
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

  # Principal Coordinates Analysis
  pmeso_pcoa <- cmdscale(sampleREL.dist,k=3,eig=TRUE,add=FALSE) 
    # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
    # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(pmeso_pcoa$eig[1]/sum(pmeso_pcoa$eig)*100,1) 
  explainvar2 <- round(pmeso_pcoa$eig[2]/sum(pmeso_pcoa$eig)*100,1)
  explainvar3 <- round(pmeso_pcoa$eig[3]/sum(pmeso_pcoa$eig)*100,1)
  
  pcoap <- merge(as.data.frame(pmeso_pcoa$points),design,by=0,all.x=T)# [,-1]
  rownames(pcoap) <- rownames(pmeso_pcoa$points)
  
 # Plot Parameters
  par(mar=c(5,5,1,1))#, oma=c(3,1,1,1)+0.1 )
  x.dim <- c(min(pcoap$V1)-(max(pcoap$V1)*0.1),max(pcoap$V1)+(max(pcoap$V1)*0.1))
  y.dim <- c(min(pcoap$V2)-(max(pcoap$V2)*0.1),max(pcoap$V2)+(max(pcoap$V2)*0.1))

  # Initiate Plot
  plot(pcoap$V1, pcoap$V2, xlab=paste("PCoA Axis 1 (",explainvar1, "%)", sep="")
    , ylab=paste("PCoA Axis 2 (",explainvar2, "%)", sep=""), 
    xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",xaxt="n",
    yaxt="n", cex.lab=1.5, cex.axis=1.2)
  axis(side=1, las=1)   
  axis(side=2, las=1)    
  abline(h=0, lty="dotted")  
  abline(v=0, lty="dotted")
  box(lwd=2)
  
  # Add Points & Ellipses
  points(pcoap$V1, pcoap$V2, pch=19, cex=2.0, bg="gray", col="gray")
  #text(pcoap$V1, pcoap$V2, labels=pcoap$Treatment, pos=3)   
  ordiellipse(cbind(pcoap$V1, pcoap$V2), pcoap$Treatment, kind="se", conf=0.95,
    lwd=2, draw = "polygon", col="gray", border = "black", label=TRUE, cex=2) 


  # Save Plots
  dev.copy2pdf(file=paste("../graphs/",plot.title,".pdf",sep=""))
  dev.copy(png, file=paste("../graphs/",plot.title,".png",sep=""), width=72*(7*4), 
    height=72*(8*4), res=72*4)
  dev.off()
  
  # Adonis (PERMANOVA)
  # Adonis runs a PERMANOVA (Created by Marti J. Anderson) 
  # this is very similar to ANOVA but for multivariate data
  # You can make very complex experimental designs with it
  # The default distance measure is bray-curtis, but other measures 
  # (Chao, Jaccard, Euclidean) can be used when specified  
  Adonis <- adonis(sampleREL.dist ~ design$Treatment, method="bray", 
    permutations=1000)
    return(Adonis)
  }
  
pmeso.pcoa.input <- function(shared = " ", design = " ", plot.title = "test",
  level = "0.03"){

  # Setup Work Environment
  source("../functions/DiversityFunctions.r")
  require("vegan")
  se <- function(x){sd(x)/sqrt(length(x))}

  # Import Site by OTU Matrix
  pmeso_data <- t(read.otu(shared, level))
  design <- read.delim(design, header=T, row.names=1)

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

  # Principal Coordinates Analysis
  pmeso_pcoa <- cmdscale(sampleREL.dist,k=3,eig=TRUE,add=FALSE)
    # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
    # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(pmeso_pcoa$eig[1]/sum(pmeso_pcoa$eig)*100,1)
  explainvar2 <- round(pmeso_pcoa$eig[2]/sum(pmeso_pcoa$eig)*100,1)
  explainvar3 <- round(pmeso_pcoa$eig[3]/sum(pmeso_pcoa$eig)*100,1)

  pcoap <- merge(as.data.frame(pmeso_pcoa$points),design,by=0,all.x=T)# [,-1]
  rownames(pcoap) <- rownames(pmeso_pcoa$points)
  pcoa <- list("pcoap" = pcoap, "var1" = explainvar1, "var2" = explainvar2,
    "var3" = explainvar3)
  return(pcoa)
  }

