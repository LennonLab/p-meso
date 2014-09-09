################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: Micro Graphs                                    #
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
bio.data <- read.delim("../data/p-meso_bio.txt", header=T, row.names=1)
design <- read.delim("../data/design.txt")

bio.data$Treatment <- factor(design$Treatment, 
  levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

################################################################################
# Plots 
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=6, height=8)
par(mfrow=c(2,1), mar=c(0.25,4.5,0.25,0.25), oma=c(3,0.5,0.5,0.5)+0.1 ) 

# Bacterial Production
bp.means <- tapply(bio.data$BP.uMC, bio.data$Treatment, mean)
bp.se <- tapply(bio.data$BP.uMC, bio.data$Treatment, se)
plot(bp.means, ylim=c(0,0.35), xlim=c(0.5, 6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab=expression(paste("BP (µM C hr"^" -1", ")", sep="")), 
  cex.lab=1.4, cex.axis=1.25)
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


# Bacterial Growth Efficiency 
bge.means <- tapply(bio.data$BGE, bio.data$Treatment, mean)
bge.se <- tapply(bio.data$BGE, bio.data$Treatment, se)
plot(bge.means, ylim=c(0.0,0.45), xlim=c(0.5,6.5), pch=15, cex=2, las=1, 
  xaxt="n", ylab="BGE (%)", cex.lab=1.4, cex.axis=1.25)
arrows(x0 = 1:6, y0 = bge.means, y1 = bge.means - bge.se, angle = 90, 
  length=0.05, lwd = 2)
arrows(x0 = 1:6, y0 = bge.means, y1 = bge.means + bge.se, angle = 90, 
  length=0.05, lwd = 2)
box(lwd=2)
axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'PA', 'Mix'), 
  at=1:6, lwd.ticks=2, cex.axis=1.2)
axis(side = 2, labels=F, lwd.ticks=2)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)
dev.copy2pdf(file=paste("./bio_plot.pdf",sep=""))
dev.copy(png, file=paste("./bio_plot.png",sep=""), bg = "white", width=72*(5*4), 
    height=72*(6*4), res=72*4)
dev.off()
windows.options(reset=TRUE)



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
