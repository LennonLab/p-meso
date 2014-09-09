################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment:  Microbial Community Diversity Graphs        #
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
setwd("~/GitHub/p-meso/graphs/")
source("../functions/p-meso.diversity.R")

# Heterotrophs
heteros <- pmeso.diversity.calc(shared = "../mothur/output/tanks.bac.final.shared",
  design = "../data/design.txt", level = "0.03", size = 10000)

# Phototrophs
photos <- pmeso.diversity.calc(shared = "../mothur/output/tanks.photo.final.shared",
  design = "../data/design.txt", level = "0.03", size = 10000)

# Eukaryotes
euks <- pmeso.diversity.calc(shared = "../mothur/output/tanks.euk.final.shared",
  design = "../data/design.txt", level = "unique", size = 100)

################################################################################
# Plots
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=10, height=6)
par(mfrow=c(2,3), mar=c(0.5,3,3,1), oma=c(1,4,1,1)+0.1 )

# Richness Plot Heteros
  rich.means <- tapply(heteros$Richness$mean, heteros$Richness$Design, mean)
  rich.se <- tapply(heteros$Richness$mean, heteros$Richness$Design, se)
  dim.r <- round(c(min(rich.means)-(max(rich.se)*1.1),
    max(rich.means)+(max(rich.se)*1.1)),2)
  plot(rich.means, ylim=c(180,300), xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  arrows(x0 = 1:6, y0 = rich.means, y1 = rich.means - rich.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = rich.means, y1 = rich.means + rich.se, angle = 90,
    length=0.05, lwd = 2)
  mtext("Bacteria", side = 3, outer=FALSE, cex=1.5, line = 1)
  mtext("Taxa Richness", side = 2, outer = TRUE, cex = 1.5, line = 1.5, adj=0.80)
  box(lwd=2)
  axis(side = 1, labels=F, lwd.ticks=2)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25, at=c(180, 205, 230, 255, 280))
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2,  at=c(180, 205, 230, 255, 280))
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2,  at=c(180, 205, 230, 255, 280))
  
# Richness Plot Photos
  rich.means <- tapply(photos$Richness$mean, photos$Richness$Design, mean)
  rich.se <- tapply(photos$Richness$mean, photos$Richness$Design, se)
  dim.r <- round(c(min(rich.means)-(max(rich.se)*1.1),
    max(rich.means)+(max(rich.se)*1.1)),2)
  plot(rich.means, ylim=c(7.5, 25), xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  arrows(x0 = 1:6, y0 = rich.means, y1 = rich.means - rich.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = rich.means, y1 = rich.means + rich.se, angle = 90,
    length=0.05, lwd = 2)
  mtext("Cyanobacteria", side = 3, outer=FALSE, cex=1.5, line = 1)
  box(lwd=2)
  axis(side = 1, labels=F, lwd.ticks=2)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)
  
# Richness Plot Euks
  rich.means <- tapply(euks$Richness$mean, euks$Richness$Design, mean)
  rich.se <- tapply(euks$Richness$mean, euks$Richness$Design, se)
  dim.r <- round(c(min(rich.means)-(max(rich.se)*1.1),
    max(rich.means)+(max(rich.se)*1.1)),2)
  plot(rich.means, ylim=dim.r, xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  arrows(x0 = 1:6, y0 = rich.means, y1 = rich.means - rich.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = rich.means, y1 = rich.means + rich.se, angle = 90,
    length=0.05, lwd = 2)
  mtext("Algae", side = 3, outer=FALSE, cex=1.5, line = 1)
  box(lwd=2)
  axis(side = 1, labels=F, lwd.ticks=2)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25, at=c(4, 8, 12, 16))
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2, at=c(4, 8, 12, 16))
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2, at=c(4, 8, 12, 16))

par(mar=c(1.5,3,0,1))
# Evenness Plot Heteros
  even.means <- tapply(heteros$Evenness$mean, heteros$Evenness$Design, mean)
  even.se <- tapply(heteros$Evenness$mean, heteros$Evenness$Design, se)
  dim.e <- round(c(min(even.means)-(max(even.se)*1.1),
    max(even.means)+(max(even.se)*1.1)),2)
  plot(even.means, ylim=dim.e, xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  arrows(x0 = 1:6, y0 = even.means, y1 = even.means - even.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = even.means, y1 = even.means + even.se, angle = 90,
    length=0.05, lwd = 2)
  mtext("Taxa Evenness", side = 2, outer = TRUE, cex = 1.5, line = 1.5, adj=0.20)
  box(lwd=2)
  axis(side = 1, labels=c('Ctrl', 'PO4', 'AEP', 'ATP', 'PA', 'Mix'), at=1:6,
    lwd.ticks=2, cex.axis=1.25)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25, at=c(0.50, 0.54, 0.58, 0.62, 0.66))
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2, at=c(0.50, 0.54, 0.58, 0.62, 0.66))
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2, at=c(0.50, 0.54, 0.58, 0.62, 0.66))
  
# Evenness Plot Photos
  even.means <- tapply(photos$Evenness$mean, photos$Evenness$Design, mean)
  even.se <- tapply(photos$Evenness$mean, photos$Evenness$Design, se)
  dim.e <- round(c(min(even.means)-(max(even.se)*1.1),
    max(even.means)+(max(even.se)*1.1)),2)
  plot(even.means, ylim=dim.e, xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  arrows(x0 = 1:6, y0 = even.means, y1 = even.means - even.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = even.means, y1 = even.means + even.se, angle = 90,
    length=0.05, lwd = 2)
  box(lwd=2)
  axis(side = 1, labels=c('Ctrl', 'PO4', 'AEP', 'ATP', 'PA', 'Mix'), at=1:6,
    lwd.ticks=2, cex.axis=1.25)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)
  
# Evenness Plot Euks
  even.means <- tapply(euks$Evenness$mean, euks$Evenness$Design, mean)
  even.se <- tapply(euks$Evenness$mean, euks$Evenness$Design, se)
  dim.e <- round(c(min(even.means)-(max(even.se)*1.1),
    max(even.means)+(max(even.se)*1.1)),2)
  plot(even.means, ylim=dim.e, xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  arrows(x0 = 1:6, y0 = even.means, y1 = even.means - even.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = even.means, y1 = even.means + even.se, angle = 90,
    length=0.05, lwd = 2)
  box(lwd=2)
  axis(side = 1, labels=c('Ctrl', 'PO4', 'AEP', 'ATP', 'PA', 'Mix'), at=1:6,
    lwd.ticks=2, cex.axis=1.25)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

  # Save Plots
  dev.copy2pdf(file=paste("../graphs/diversity.pdf",sep=""))
  dev.copy(png, file=paste("../graphs/diversity.png",sep=""),
    bg="white", width=72*(4*4), height=72*(8*4), res=72*4)
  dev.off()
  windows.options(reset=TRUE)
