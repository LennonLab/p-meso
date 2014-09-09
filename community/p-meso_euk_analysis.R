################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment: Eukaryote PcoA and PERMANOVA                 #
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
tanks_0.03 <- t(read.otu("../mothur/output/tanks.euks.shared", "0.03"))
tanks_0.05 <- t(read.otu("../mothur.output/tanks.euks.shared", "0.05"))

# Make Presence Absence Matrices
Tanks_0.03PA <- (Tanks_0.03>0)*1
Tanks_0.05PA <- (Tanks_0.05>0)*1
# Make Relative Abundence Matrices
Tanks_0.03REL <- Tanks_0.03
Tanks_0.05REL <- Tanks_0.05
for(i in 1:24){
	Tanks_0.03REL[,i]<-Tanks_0.03[,i]/sum(Tanks_0.03[,i])
	Tanks_0.05REL[,i]<-Tanks_0.05[,i]/sum(Tanks_0.05[,i])
	}

# I am doing a few things here. I am using decostand (in vegan) to log transform the abundance data. I am using a rel abund matrix so it first divides by the lowest value. I am using this transformation to help with the bias against low abundance (rare species). See Anderson et al 2006 or Legendre & Legendre 2012 p 327
Tanks_0.03REL.dist <- vegdist(decostand(t(Tanks_0.03REL), method="log"),method="bray")
Tanks_0.03REL.dist

# visualization of matrix
Tanks_0.03REL.pcoa <- cmdscale(Tanks_0.03REL.dist, k=2)
Tanks_0.03REL.pcoa
Tanks_Pcoa <- as.data.frame(Tanks_0.03REL.pcoa)
Tanks_Pcoa <- merge(Tanks_Pcoa, design, by=0)
colnames(Tanks_Pcoa) <- c("Tank", "Axis 1", "Axis 2", "Treatment")
plot(Tanks_Pcoa[,2], Tanks_Pcoa[,3], type='n')
text(Tanks_Pcoa[,2], Tanks_Pcoa[,3], labels=Tanks_Pcoa$Treatment)

par(mar=c(5,5,1,1))#, oma=c(3,1,1,1)+0.1 )
plot(Tanks_0.03REL.pcoa,  bg="black", col="black", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.2, xlab="PCoA 1", ylab="PCoA 2", xlim=c(-0.27,0.3), ylim=c(-0.21,0.33),las=1)
text(Tanks_0.03REL.pcoa, labels=Tanks_Pcoa$Treatment, pos=3)
ordiellipse(Tanks_0.03REL.pcoa, Tanks_Pcoa$Treatment, kind="se", conf=0.95, lwd=2, draw = "polygon", col="skyblue", border = "blue", label=TRUE)
box(lwd=2)

# Black and White
par(mar=c(5,5,1,1))#, oma=c(3,1,1,1)+0.1 )
plot(Tanks_0.03REL.pcoa,  bg="black", col="black", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.2, xlab="PCoA 1", ylab="PCoA 2", xlim=c(-0.27,0.3), ylim=c(-0.21,0.33),las=1)
text(Tanks_0.03REL.pcoa, labels=Tanks_Pcoa$Treatment, pos=3)
ordiellipse(Tanks_0.03REL.pcoa, Tanks_Pcoa$Treatment, kind="se", conf=0.95, lwd=2, draw = "polygon", col="gray", border = "black", label=TRUE)
box(lwd=2)


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



# Alpha Diversity with Resampling
Rich <- round(richness.iter(shared = "tanks.euks.shared", cutoff = "0.03", size = 10000, iters = 100), 3)
Even <- round(evenness.iter(shared = "tanks.euks.shared", cutoff = "0.03", size = 10000, iters = 100), 3)
Shan <- round(diversity.iter(shared = "tanks.euks.shared", index = "shannon", cutoff = "0.03", size = 10000, iters = 100), 3)

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
plot(Rich_data$Design, Rich_data$mean, type='n', ylim=c(22, 38), yaxt="n", cex=1.5, cex.lab=2, cex.axis=1.2, las=1, border="white", xaxt="n")
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
plot(Even_data$Design, Even_data$mean, type='n', ylim=c(0.15,0.45), yaxt="n", cex=1.5, cex.lab=2, cex.axis=1.2, las=1, border="white", xaxt="n")
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
plot(Shan_data$Design, Shan_data$mean, type='n', ylim=c(0.4,1.5), yaxt="n", cex=1.5, cex.lab=2, cex.axis=1.2, las=1, border="white", xaxt="n")
points(Shan.means, pch=15, cex=2)
arrows(x0 = 1:6, y0 = Shan.means, y1 = Shan.means - Shan.se, angle = 90, length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = Shan.means, y1 = Shan.means + Shan.se, angle = 90, length=0.05, lwd = 2)
title(ylab = "Shannon Diversity", cex.lab = 2, line = 3.5)
box(lwd=2)
axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'), at=1:6, lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
axis(side = 3, labels=F, lwd.ticks=2)
axis(side = 4, labels=F, lwd.ticks=2)