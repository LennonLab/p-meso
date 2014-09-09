################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: Foodweb & Eco Graph                             #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/06/26                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())  
setwd("~./GitHub/p-meso/graphs")

# Import Data
zoop.data <- read.delim("../data/zoop.txt", header=T, row.names=1)
eco.data <- read.delim("../data/p-meso_eco.txt", header=T, row.names=1)
design <- read.delim("../data/design.txt")

design$Treatment <- factor(design$Treatment, 
  levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix')) 
  zoop.data$Treatment <- design$Treatment 
  eco.data$Treatment <- design$Treatment 
  
# Calculations
zoop.data$Cladocera_S <- round(((zoop.data$Cladocera * (zoop.data$Vol / zoop.data$Vol.subSamp)) / zoop.data$Tows),1)
zoop.data$Copepoda_S <- round(((zoop.data$Copepoda * (zoop.data$Vol / zoop.data$Vol.subSamp)) / zoop.data$Tows),1)
zoop.data$Total_S <- zoop.data$Cladocera_S + zoop.data$Copepoda_S
zoop.data$CladCope <- zoop.data$Cladocera_S / zoop.data$Copepoda_S

################################################################################
# Plots 
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=9, height=6)
par(mfrow=c(2,2), mar=c(0.25,4.5,0.25,1), oma=c(3,1,1,1)+0.1 )

# Chl A
ch.means <- tapply(eco.data$ChlA2, eco.data$Treatment, mean)
ch.se <- tapply(eco.data$ChlA2, eco.data$Treatment, se)
plot(ch.means, ylim=c(0,6), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  axes=FALSE, ylab=expression(paste("Chl ", italic("a")," (µg l"^" -1",")",
  sep="")), cex.lab=1.5, cex.axis=1.25)
arrows(x0 = 1:6, y0 = ch.means, y1 = ch.means - ch.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = ch.means, y1 = ch.means + ch.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
text(0.6, 5.7, "A", bty="n", cex=1.5)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25, at=c(0,2,4,6))
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2, at=c(0,2,4,6))
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2, at=c(0,2,4,6))

# Respiration
rp.means <- tapply(eco.data$Resp_uM, eco.data$Treatment, mean, na.rm=TRUE)
rp.se <- tapply(eco.data$Resp_uM, eco.data$Treatment, se)
plot(rp.means, ylim=c(0,50), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("R (µmol C l"^" -1", " hr"^" -1",")",
  sep="")), cex.lab=1.5, cex.axis=1.25)
arrows(x0 = 1:6, y0 = rp.means, y1 = rp.means - rp.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = rp.means, y1 = rp.means + rp.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
text(0.6, 47.5, "C", bty="n", cex=1.5)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Clad:Cope
ratio.means <- tapply(zoop.data$CladCope, zoop.data$Treatment, mean)
ratio.se <- tapply(zoop.data$CladCope, zoop.data$Treatment, se)
plot(ratio.means, ylim=c(0,3), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  axes=FALSE, ylab="Clad:Cope", cex.lab=1.5, cex.axis=1.25)
arrows(x0 = 1:6, y0 = ratio.means, y1 = ratio.means - ratio.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = ratio.means, y1 = ratio.means + ratio.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
text(0.6, 2.85, "B", bty="n", cex=1.5)
axis(side = 1, labels=c('Ctrl', 'PO4', 'AEP', 'ATP', 'PA', 'Mix'), 
  at=1:6, lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25, at=c(0,1,2,3))
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2, at=c(0,1,2,3))
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2, at=c(0,1,2,3))

# NPP
npp.means <- tapply(eco.data$NPP_uM, eco.data$Treatment, mean)
npp.se <- tapply(eco.data$NPP_uM, eco.data$Treatment, se)
plot(npp.means, ylim=c(25,65), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("NPP (µmol C l"^" -1", " hr"^" -1",")"
  ,sep="")), cex.lab=1.5, cex.axis=1.25)
arrows(x0 = 1:6, y0 = npp.means, y1 = npp.means - npp.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = npp.means, y1 = npp.means + npp.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
text(0.6, 63, "D", bty="n", cex=1.5)
axis(side = 1, labels=c('Ctrl', 'PO4', 'AEP', 'ATP', 'PA', 'Mix'), 
  at=1:6, lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Save Plots
dev.copy2pdf(file=paste("./foodweb-eco_plot.pdf",sep=""))
dev.copy(png, file=paste("./foodweb-eco_plot.png",sep=""), bg="white", width=72*(6*4), 
    height=72*(4*4), res=72*4)
dev.off()
windows.options(reset=TRUE)

