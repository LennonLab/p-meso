################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment: Heterotroph PcoA and PERMANOVA               #
#                                                                              #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/24                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/p-meso/community/")
source("../functions/DiversityFunctions.r")
require("vegan")
require("ecodist")
se <- function(x){sd(x)/sqrt(length(x))}

# Import Data
design <- read.delim("../data/design.txt", row.names=1)
tanks_0.03 <- t(read.otu("../mothur/output/tanks.bac.final.shared", "0.03"))
tanks_0.05 <- t(read.otu("../mothur/output/tanks.bac.final.shared", "0.05"))

# Make Presence Absence Matrices
tanks_0.03PA <- (tanks_0.03>0)*1
tanks_0.05PA <- (tanks_0.05>0)*1

# Make Relative Abundence Matrices
tanks_0.03REL <- tanks_0.03
  for(i in 1:ncol(tanks_0.03REL)){
    tanks_0.03REL[,i] = tanks_0.03[,i]/sum(tanks_0.03[,i])
    } 
tanks_0.05REL <- tanks_0.05
  for(i in 1:ncol(tanks_0.05REL)){
    tanks_0.05REL[,i] = tanks_0.05[,i]/sum(tanks_0.05[,i])
    }  
    
# Create Distance Matrix
tanks_PA.dist <- vegdist(decostand(t(tanks_0.03PA),method="log"),
  method="bray")
tanks_REL.dist <- vegdist(decostand(t(tanks_0.03REL),method="log"),
  method="bray")
  # I am using decostand (in vegan) to log transform the abundance data. 
  # I am using a rel abund matrix so it first divides by the lowest value. 
  # I am using this transformation to help with the bias against low abundance 
  # See Anderson et al 2006 or Legendre & Legendre 2012 p 327

# Principal Coordinates Analysis
tanks_pcoa <- cmdscale(tanks_REL.dist,k=3,eig=TRUE,add=FALSE) 
  # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
  # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

# Percent Variance Explained Using PCoA (Axis 1,2,3)
explainvar1 <- round(tanks_pcoa$eig[1]/sum(tanks_pcoa$eig)*100,2) 
explainvar2 <- round(tanks_pcoa$eig[2]/sum(tanks_pcoa$eig)*100,2)
explainvar3 <- round(tanks_pcoa$eig[3]/sum(tanks_pcoa$eig)*100,2)
  
pcoap <- merge(as.data.frame(tanks_pcoa$points),design,by=0,all.x=T)[,-1]
rownames(pcoap) <- rownames(tanks_pcoa$points)


# Plot Parameters
par(mar=c(5,5,1,1))#, oma=c(3,1,1,1)+0.1 )
x.dim <- c(min(pcoap$V2)-min(pcoap$V2)*0.2,max(pcoap$V2)+max(pcoap$V2)*0.2)
y.dim <- c(min(pcoap$V1)-min(pcoap$V1)*0.2,max(pcoap$V1)+max(pcoap$V1)*0.2)

# Initiate Plot
plot(Tanks_0.03REL.pcoa, xlab=paste("PCoA Axis 2 (",explainvar2, "%)", sep=""),
  ylab=paste("PCoA Axis 1 (",explainvar1, "%)", sep=""), xlim=x.dim, ylim=y.dim
  bg="black", col="black", pch=19, cex=2.0, cex.lab=1.5, cex.axis=1.2, las=1)
axis(side=1, las=1)   
axis(side=2, las=1)    
abline(h=0, lty="dotted")  
abline(v=0, lty="dotted")
text(Tanks_0.03REL.pcoa, labels=Tanks_Pcoa$Treatment, pos=3)
ordiellipse(Tanks_0.03REL.pcoa, Tanks_Pcoa$Treatment, kind="se", conf=0.95, 
  lwd=2, draw = "polygon", col="skyblue", border = "blue", label=TRUE)
box(lwd=2)
dev.copy2pdf(file=paste("./plots/",plot.title,".pdf",sep=""))
dev.copy(png, file=paste("./plots/",plot.title,".png",sep=""), width=72*(7*4), 
  height=72*(8*4), res=72*4)
dev.off()



# Black and White
par(mar=c(5,5,1,1))#, oma=c(3,1,1,1)+0.1 )
plot(Tanks_0.03REL.pcoa,  bg="black", col="black", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.2, xlab="PCoA 1", ylab="PCoA 2", xlim=c(-0.25,0.2), ylim=c(-0.35,0.2),las=1)
text(Tanks_0.03REL.pcoa, labels=Tanks_Pcoa$Treatment, pos=3)
ordiellipse(Tanks_0.03REL.pcoa, Tanks_Pcoa$Treatment, kind="se", conf=0.95, lwd=2, draw = "polygon", col="gray", border = "black", label=TRUE)
box(lwd=2)


  
  # Initiate Plot
  plot(pcoap$V2, pcoap$V1, xlab=paste("PCoA Axis 2 (",explainvar2, "%)", sep="")
    , ylab=paste("PCoA Axis 1 (",explainvar1, "%)", sep=""), 
    xlim=x.dim,ylim= y.dim, pch=16, cex=2.0, type="n",xaxt="n",
    yaxt="n", cex.lab=1.5, cex.axis=1.2)  
  axis(side=1, las=1)   
  axis(side=2, las=1)    
  abline(h=0, lty="dotted")  
  abline(v=0, lty="dotted")
  mol.shape <- rep(NA, dim(pcoap)[1])
    for (i in 1:length(mol.shape)){
      if (pcoap$molecule[i] == "DNA"){mol.shape[i] = 21}
      else {mol.shape[i] = 22}
      }
  slope.color <- rep(NA, dim(pcoap)[1])
    for (i in 1:length(slope.color)){
      if (pcoap$slope[i] == levels(pcoap$slope)[1]) {slope.color[i] = "brown2"}
      else {slope.color[i] = "green3"}
      } 
  points(pcoap$V2, pcoap$V1, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)   
  ordiellipse(cbind(pcoap$V2, pcoap$V1), pcoap$site, kind="sd", conf=0.95,
    lwd=2, lty=3, draw = "lines", col = "black", label=TRUE)
  # legend("topleft", c(paste("All; ",levels(pcoap$slope)[1]," Slope", sep=""), 
  #   paste("All; ",levels(pcoap$slope)[2]," Slope", sep=""), 
  #   paste("Active; ",levels(pcoap$slope)[1]," Slope", sep=""),
  #   paste("Active; ",levels(pcoap$slope)[2]," Slope", sep="")), 
  #   pt.lwd=2, col="black", pt.bg=c("brown3", "green3", "brown3", 
  #   "green3"), pch=c(21,21,22,22), bty='o', box.lty=0, bg="white", cex=1.5)
  box(lwd=2)
  par(mar=c(0, 3, 0, 0))
  plot.new()
  legend("center", c(paste("All; ",levels(pcoap$slope)[1]," Slope", sep=""), 
    paste("All; ",levels(pcoap$slope)[2]," Slope", sep=""), 
    paste("Active; ",levels(pcoap$slope)[1]," Slope", sep=""),
    paste("Active; ",levels(pcoap$slope)[2]," Slope", sep="")), 
    pt.lwd=2, col="black", pt.bg=c("brown2", "green3", "brown2", 
    "green3"), pch=c(21,21,22,22), bty='n', ncol=2, cex=1.5, pt.cex=2)
  dev.copy2pdf(file=paste("./plots/",plot.title,".pdf",sep=""))
  dev.copy(png, file=paste("./plots/",plot.title,".png",sep=""), width=72*(7*4), 
    height=72*(8*4), res=72*4)
  dev.off()
  
  # Adonis (PERMANOVA)
  # Adonis runs a PERMANOVA (Created by Marti J. Anderson) 
  # this is very similar to ANOVA but for multivariate data
  # You can make very complex experimental designs with it
  # The default distance measure is bray-curtis, but other measures 
  # (Chao, Jaccard, Euclidean) can be used when specified  
#  Adonis <- adonis(sampleREL.dist ~ design$molecule*design$slope, method="bray", 
#    permutations=1000)
#    return(Adonis)
  }
  



# adonis runs a permanova (Created by Marti J. Anderson) this is very similar to ANOVA but for multivariate data. You can make very complex
# experimental designs with it.
# The default distance measure is bray-curtis, but other measures (Chao, Jaccard, Euclidean) can be used when specified
adonis(Tanks_0.03REL.dist ~ design$Treatment , method="bray", permutations=1000)



# I am using metaMDS (vegan) to do NonMetric Multi Demensional Scalling. 
Tanks_0.03REL.nmds <- metaMDS(Tanks_0.03REL.dist, k=2)
par(mar=c(5,5,1,1))
plot(Tanks_0.03REL.nmds, display="sites", xlab="NMDS 1", ylab="NMDS 2",las=1, cex.lab=1.5, cex.axis=1.2)
text(Tanks_0.03REL.nmds, labels=design$Treatment, pos=3)
ordiellipse(Tanks_0.03REL.nmds, design$Treatment, kind="se", conf=0.95, lwd=2, draw = "polygon", col="skyblue", border = "blue", label=TRUE)
box(lwd=2)


# Alpha Diversity Measures
Tanks_Rich <- richness.iter("Tanks.shared", "0.03", size = 10000, iters = 100)

# Alpha Diversity with Resampling
Rich <- round(richness.iter(shared = "Tanks.shared", cutoff = "0.03", size = 10000, iters = 100), 3)
Even <- round(evenness.iter(shared = "Tanks.shared", cutoff = "0.03", size = 10000, iters = 100), 3)
Shan <- round(diversity.iter(shared = "Tanks.shared", index = "shannon", cutoff = "0.03", size = 10000, iters = 100), 3)

# Richness Summary
Rich_data <- merge(design, Rich, by="row.names")
Rich_data$mean <- round(apply(Rich_data[3:26], 1, mean),3)
Rich_data$se <- round(apply(Rich_data[3:26], 1, se), 3)
Rich_data$Design <- design$Treatment
Rich_data$Design <- factor(Rich_data$Design, levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

# Evenness Summary
Even_data <- merge(design, Even, by="row.names")
Even_data$mean <- round(apply(Even_data[3:26], 1, mean),3)
Even_data$se <- round(apply(Even_data[3:26], 1, se),3)
Even_data$Design <- design$Treatment
Even_data$Design <- factor(Even_data$Design, levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

# Shannon Diversity Summary
Shan_data <- merge(design, Shan, by="row.names")
Shan_data$mean <- round(apply(Shan_data[3:26], 1, mean),3)
Shan_data$se <- round(apply(Shan_data[3:26], 1, se),3)
Shan_data$Design <- design$Treatment
Shan_data$Design <- factor(Shan_data$Design, levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

#par(mar=c(3,5,1,1))
windows.options(width=6, height=12)
par(mfrow=c(3,1), mar=c(0.25,5,0.25,1), oma=c(3,1,1,1)+0.1 )

# Richness Plot A
Rich.means <- tapply(Rich_data$mean, Rich_data$Design, mean)
Rich.se <- tapply(Rich_data$mean, Rich_data$Design, se)
plot(Rich_data$Design, Rich_data$mean, type='n', ylim=c(180, 300), yaxt="n", cex=1.5, cex.lab=2, cex.axis=1.2, las=1, border="white", xaxt="n")
points(Rich.means, pch=15, cex=2)
arrows(x0 = 1:6, y0 = Rich.means, y1 = Rich.means - Rich.se, angle = 90, length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = Rich.means, y1 = Rich.means + Rich.se, angle = 90, length=0.05, lwd = 2)
title(ylab = "Taxa Richness", cex.lab = 2, line = 3.5)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
axis(side = 3, labels=F, lwd.ticks=2)
axis(side = 4, labels=F, lwd.ticks=2)

# Evenness Plot A
Even.means <- tapply(Even_data$mean, Even_data$Design, mean)
Even.se <- tapply(Even_data$mean, Even_data$Design, se)
plot(Even_data$Design, Even_data$mean, type='n', ylim=c(0.53,0.68), yaxt="n", cex=1.5, cex.lab=2, cex.axis=1.2, las=1, border="white", xaxt="n")
points(Even.means, pch=15, cex=2)
arrows(x0 = 1:6, y0 = Even.means, y1 = Even.means - Even.se, angle = 90, length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = Even.means, y1 = Even.means + Even.se, angle = 90, length=0.05, lwd = 2)
title(ylab = "Taxa Evenness", cex.lab = 2, line = 3.5)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
axis(side = 3, labels=F, lwd.ticks=2)
axis(side = 4, labels=F, lwd.ticks=2)

# Shannon Plot A
Shan.means <- tapply(Shan_data$mean, Shan_data$Design, mean)
Shan.se <- tapply(Shan_data$mean, Shan_data$Design, se)
plot(Shan_data$Design, Shan_data$mean, type='n', ylim=c(2.8,3.8), yaxt="n", cex=1.5, cex.lab=2, cex.axis=1.2, las=1, border="white", xaxt="n")
points(Shan.means, pch=15, cex=2)
arrows(x0 = 1:6, y0 = Shan.means, y1 = Shan.means - Shan.se, angle = 90, length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = Shan.means, y1 = Shan.means + Shan.se, angle = 90, length=0.05, lwd = 2)
title(ylab = "Shannon Diversity", cex.lab = 2, line = 3.5)
box(lwd=2)
axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'), at=1:6, lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
axis(side = 3, labels=F, lwd.ticks=2)
axis(side = 4, labels=F, lwd.ticks=2)

windows.options(reset=TRUE)