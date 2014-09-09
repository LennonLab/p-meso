################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: Eco Graphs                                     #
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
eco.data <- read.delim("../data/p-meso_eco.txt", header=T, row.names=1)
design <- read.delim("../data/design.txt")

eco.data$Treatment <- factor(design$Treatment, 
  levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

################################################################################
# Plots 
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=6, height=8)
par(mfrow=c(2,2), mar=c(0.25,4.5,0.25,0.25), oma=c(3,0.5,0.5,0.5)+0.1 ) 

# Chl A
ch.means <- tapply(eco.data$ChlA2, eco.data$Treatment, mean)
ch.se <- tapply(eco.data$ChlA2, eco.data$Treatment, se)
plot(ch.means, ylim=c(0,6), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("Chl ", italic("a")," (µg L"^"-1",")",
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

# NPP
npp.means <- tapply(eco.data$NPP_uM, eco.data$Treatment, mean)
npp.se <- tapply(eco.data$NPP_uM, eco.data$Treatment, se)
plot(npp.means, ylim=c(25,65), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("NPP (µM O"["2"]," hr"^"-1",
  ")",sep="")), cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = npp.means, y1 = npp.means - npp.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = npp.means, y1 = npp.means + npp.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'PA', 'Mix'), 
  at=1:6, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)
dev.copy2pdf(file=paste("./eco_plot.pdf",sep=""))
dev.copy(png, file=paste("./eco_plot.png",sep=""), bg = "white", width=72*(5*4), 
    height=72*(6*4), res=72*4)
dev.off()
windows.options(reset=TRUE)



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
