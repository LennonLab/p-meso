################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: 3 Panel Graph                                   #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/17                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())  
setwd("~./GitHub/p-meso/graphs")

# Import Data
chem.data <- read.delim("../data/p-meso_chem.txt", header=T, row.names=1)
bio.data <- read.delim("../data/p-meso_bio.txt", header=T, row.names=1)
eco.data <- read.delim("../data/p-meso_eco.txt", header=T, row.names=1)
design <- read.delim("../data/design.txt")

design$Treatment <- factor(design$Treatment, 
  levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))
  chem.data$Treatment <- design$Treatment 
  bio.data$Treatment <- design$Treatment 
  eco.data$Treatment <- design$Treatment 

################################################################################
# Plots 
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=18, height=12)
par(mfrow=c(3,3), mar=c(0.25,4.5,0.25,1), oma=c(3,1,1,1)+0.1 )

# Total Phosphorous
tp.means <- tapply(chem.data$TP.uM, chem.data$Treatment, mean)
tp.se <- tapply(chem.data$TP.uM, chem.data$Treatment, se)
plot(tp.means, ylim=c(0,0.6), xlim = c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Total Phosphorus (µM L"^"-1",")", sep="")), 
  cex.lab=1.4, cex.axis=1.2)
arrows(x0 = 1:6, y0 = tp.means, y1 = tp.means - tp.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = tp.means, y1 = tp.means + tp.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Bacterial Production
bp.means <- tapply(bio.data$BP.uMC, bio.data$Treatment, mean)
bp.se <- tapply(bio.data$BP.uMC, bio.data$Treatment, se)
plot(bp.means, ylim=c(0,0.35), xlim=c(0.5, 6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Bac. Production (µM C hr"^" -1", ")", sep="")), 
  cex.lab=1.5, cex.axis=1.25)
arrows(x0 = 1:6, y0 = bp.means, y1 = bp.means - bp.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = bp.means, y1 = bp.means + bp.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Chl A
ch.means <- tapply(eco.data$ChlA2, eco.data$Treatment, mean)
ch.se <- tapply(eco.data$ChlA2, eco.data$Treatment, se)
plot(ch.means, ylim=c(0,6), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Chlorophyll ", italic("a")," (µg L"^"-1",")",
  sep="")), cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = ch.means, y1 = ch.means - ch.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = ch.means, y1 = ch.means + ch.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Total  Nitrogen
tn.means <- tapply(chem.data$TN.uM, chem.data$Treatment, mean, na.rm=TRUE)
tn.se <- tapply(chem.data$TN.uM, chem.data$Treatment, se)
plot(tn.means, ylim=c(3,10), xlim = c(0.5,6.5), pch=15, cex=2, las=1,  
  xaxt="n", ylab=expression(paste("Total Dissolved Nitrogen (µM L"^"-1",")", sep="")), 
  cex.lab=1.4, cex.axis=1.2)
arrows(x0 = 1:6, y0 = tn.means, y1 = tn.means - tn.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = tn.means, y1 = tn.means + tn.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Bacterial Respiration
br.means <- tapply(bio.data$BR.uMO2, bio.data$Treatment, mean, na.rm=TRUE)
br.se <- tapply(bio.data$BR.uMO2, bio.data$Treatment, se)
plot(br.means, ylim=c(0.3,1.3), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Bac. Respiration (µM O"["2"]," hr"^" -1",")", 
  sep="")), cex.lab=1.5, cex.axis=1.25)
arrows(x0 = 1:6, y0 = br.means, y1 = br.means - br.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = br.means, y1 = br.means + br.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Respiration
rp.means <- tapply(eco.data$Resp_uM, eco.data$Treatment, mean, na.rm=TRUE)
rp.se <- tapply(eco.data$Resp_uM, eco.data$Treatment, se)
plot(rp.means, ylim=c(0,50), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Respiration (µM O"["2"]," hr"^"-1",")",
  sep="")), cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = rp.means, y1 = rp.means - rp.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = rp.means, y1 = rp.means + rp.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Dissolved Organic Carbon
doc.means <- tapply(chem.data$DOC.uM, chem.data$Treatment, mean)
doc.se <- tapply(chem.data$DOC.uM, chem.data$Treatment, se)
plot(doc.means, ylim=c(100,250), xlim = c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Dissolved Organic Carbon (µM L"^"-1",")", 
  sep="")), cex.lab=1.4, cex.axis=1.2)
arrows(x0 = 1:6, y0 = doc.means, y1 = doc.means - doc.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = doc.means, y1 = doc.means + doc.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'), 
  at=1:6, lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Bacterial Growth Efficiency 
bge.means <- tapply(bio.data$BGE, bio.data$Treatment, mean)
bge.se <- tapply(bio.data$BGE, bio.data$Treatment, se)
plot(bge.means, ylim=c(0.0,0.45), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab="Bac. Growth Efficiency (%)", cex.lab=1.5, cex.axis=1.25)
arrows(x0 = 1:6, y0 = bge.means, y1 = bge.means - bge.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = bge.means, y1 = bge.means + bge.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'), 
  at=1:6, lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# NPP
npp.means <- tapply(eco.data$NPP_uM, eco.data$Treatment, mean)
npp.se <- tapply(eco.data$NPP_uM, eco.data$Treatment, se)
plot(npp.means, ylim=c(25,65), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Net 1° Productivity (µM O"["2"]," hr"^"-1",
  ")",sep="")), cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = npp.means, y1 = npp.means - npp.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = npp.means, y1 = npp.means + npp.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'), 
  at=1:6, lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Save Plots
dev.copy2pdf(file=paste("./multi_plot.pdf",sep=""))
dev.copy(png, file=paste("./multi_plot.png",sep=""), bg="white", width=72*(9*4), 
    height=72*(6*4), res=72*4)
dev.off()
windows.options(reset=TRUE)

