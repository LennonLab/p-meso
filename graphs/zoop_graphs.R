################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: Zooplankton Graphs                              #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/04/11                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())  
setwd("~./GitHub/p-meso/graphs")

# Import Data
zoop.data <- read.delim("../data/zoop.txt", header=T, row.names=1)
design <- read.delim("../data/design.txt")

zoop.data$Treatment <- factor(design$Treatment, 
  levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))
  
# Calculations
zoop.data$Cladocera_S <- round(((zoop.data$Cladocera * (zoop.data$Vol / zoop.data$Vol.subSamp)) / zoop.data$Tows),1)
zoop.data$Copepoda_S <- round(((zoop.data$Copepoda * (zoop.data$Vol / zoop.data$Vol.subSamp)) / zoop.data$Tows),1)
zoop.data$Total_S <- zoop.data$Cladocera_S + zoop.data$Copepoda_S
zoop.data$CladCope <- zoop.data$Cladocera_S / zoop.data$Copepoda_S

################################################################################
# Plots 
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=6, height=12)
par(mfrow=c(3,1), mar=c(0.25,4.5,0.25,1), oma=c(3,1,1,1)+0.1 ) 

# Total
total.means <- tapply(zoop.data$Total_S, zoop.data$Treatment, mean)
total.se <- tapply(zoop.data$Total_S, zoop.data$Treatment, se)
plot(total.means, ylim=c(100,1250), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Total Animals L"^"-1",
  sep="")), cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = total.means, y1 = total.means - total.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = total.means, y1 = total.means + total.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Clad:Cope
ratio.means <- tapply(zoop.data$CladCope, zoop.data$Treatment, mean)
ratio.se <- tapply(zoop.data$CladCope, zoop.data$Treatment, se)
plot(ratio.means, ylim=c(0,5), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab="Cladocera:Copepoda", cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = ratio.means, y1 = ratio.means - ratio.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = ratio.means, y1 = ratio.means + ratio.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Cladocera
clad.means <- tapply(zoop.data$Cladocera_S, zoop.data$Treatment, mean, na.rm=TRUE)
clad.se <- tapply(zoop.data$Cladocera_S, zoop.data$Treatment, se)
plot(clad.means, ylim=c(20,700), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Cladocera (Animals L"^"-1",")",
  sep="")), cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = clad.means, y1 = clad.means - clad.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = clad.means, y1 = clad.means + clad.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Copepoda
cope.means <- tapply(zoop.data$Copepoda_S, zoop.data$Treatment, mean)
cope.se <- tapply(zoop.data$Copepoda_S, zoop.data$Treatment, se)
plot(cope.means, ylim=c(50,550), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Copepoda (Animals L"^"-1",")",
  sep="")), cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = cope.means, y1 = cope.means - cope.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = cope.means, y1 = cope.means + cope.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'), 
  at=1:6, lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)
dev.copy2pdf(file=paste("./zoop_plot.pdf",sep=""))
dev.copy(png, file=paste("./zoop_plot.png",sep=""), bg = "white", width=72*(4*4), 
    height=72*(8*4), res=72*4)
dev.off()
windows.options(reset=TRUE)
