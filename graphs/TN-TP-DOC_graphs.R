################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: Chem Graphs                                     #
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
design <- read.delim("../data/design.txt")

chem.data$Treatment <- factor(design$Treatment, 
  levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

################################################################################
# Plots in mg ##################################################################
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=6, height=12)
par(mfrow=c(3,1), mar=c(0.25,4.5,0.25,1), oma=c(3,1,1,1)+0.1 )
 
# Total Phosphorous
tp.means <- tapply(chem.data$TP.ug, chem.data$Treatment, mean, na.rm=TRUE)
tp.se <- tapply(chem.data$TP.ug, chem.data$Treatment, se)
plot(tp.means, ylim=c(0,20), xlim = c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Total Phosphorus (µg L"^"-1",")", sep="")), 
  cex.lab=1.4, cex.axis=1.25)
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

# Total  Nitrogen
tn.means <- tapply(chem.data$TN.mg, chem.data$Treatment, mean, na.rm=TRUE)
tn.se <- tapply(chem.data$TN.mg, chem.data$Treatment, se)
plot(tn.means, ylim=c(0.2,0.6), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Total Nitrogen (mg L"^"-1",")", sep="")), 
  cex.lab=1.4, cex.axis=1.25)
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

# Dissolved Organic Carbon
doc.means <- tapply(chem.data$DOC.mg, chem.data$Treatment, mean)
doc.se <- tapply(chem.data$DOC.mg, chem.data$Treatment, se)
plot(doc.means, ylim=c(1.3,3.2), xlim=c(0.5,6.5), pch=15, cex=2, las=1,
   xaxt="n", ylab=expression(paste("Dissolved Organic Carbon (mg L"^"-1",")", 
   sep="")), cex.lab=1.4, cex.axis=1.25)
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
dev.copy2pdf(file=paste("./chem.g_plot.pdf",sep=""))
dev.copy(png, file=paste("./chem.g_plot.png",sep=""),bg="white", width=72*(4*4), 
    height=72*(8*4), res=72*4)
dev.off()
windows.options(reset=TRUE)

################################################################################
# Plots in µM ##################################################################
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=6, height=12)
par(mfrow=c(3,1), mar=c(0.25,4.5,0.25,1), oma=c(3,1,1,1)+0.1 ) 

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

# Total  Nitrogen
tn.means <- tapply(chem.data$TN.uM, chem.data$Treatment, mean, na.rm=TRUE)
tn.se <- tapply(chem.data$TN.uM, chem.data$Treatment, se)
plot(tn.means, ylim=c(3,10), xlim = c(0.5,6.5), pch=15, cex=2, las=1,  
  xaxt="n", ylab=expression(paste("Total Nitrogen (µM L"^"-1",")", sep="")), 
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
dev.copy2pdf(file=paste("./chem.m_plot.pdf",sep=""))
dev.copy(png, file=paste("./chem.m_plot.png",sep=""),bg="white", width=72*(4*4), 
    height=72*(8*4), res=72*4)
dev.off()
windows.options(reset=TRUE)