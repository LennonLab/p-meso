################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment:  Microbial Community PCoA Graphs             #
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
setwd("~/GitHub/p-meso/community/")
source("../functions/p-meso.pcoa.R")

# Heterotrophs
heteros <- pmeso.pcoa.input(shared = "../mothur/output/tanks.bac.final.shared",
  design = "../data/design2.txt", plot.title = "heterotroph_pcoa", level = "0.03")

# Phototrophs
photos <- pmeso.pcoa.input(shared = "../mothur/output/tanks.photo.final.shared",
  design = "../data/design2.txt", plot.title = "phototroph_pcoa", level = "0.03")

# Eukaryotes
euks <- pmeso.pcoa.input(shared = "../mothur/output/tanks.euk.final.shared",
  design = "../data/design2.txt", plot.title = "eukaryote_pcoa", level = "0.03")

################################################################################
# Plots
################################################################################
se <- function(x){x = na.omit(x);sd(x)/sqrt(length(x))}
windows.options(width=4, height=10)
par(mfrow=c(3,1), mar=c(4.5,4.5,1,1), oma=c(1,1,1,1)+0.1 )
          
# Heteros
  x.dim <- c(min(heteros$pcoap$V1)-(max(heteros$pcoap$V1)*0.1),
    max(heteros$pcoap$V1)+(max(heteros$pcoap$V1)*0.1))
  y.dim <- c(min(heteros$pcoap$V2)-(max(heteros$pcoap$V2)*0.1),
    max(heteros$pcoap$V2)+(max(heteros$pcoap$V2)*0.1))
  # Initiate Plot
  plot(heteros$pcoap$V1, heteros$pcoap$V2,
    xlab=paste("PCoA 1 (",heteros$var1, "%)", sep=""),
    ylab=paste("PCoA 2 (",heteros$var2, "%)", sep=""),
    xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",
    cex.lab=1.5, cex.axis=1.2, axes=FALSE)
  axis(side = 1, labels=T, lwd.ticks=2, at=c(-0.3,-0.15,0,0.15), cex.axis=1.2)
  axis(side = 2, labels=T, lwd.ticks=2, at=c(-0.3,-0.15,0,0.15), cex.axis=1.2)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  #text(min(heteros$pcoap$V1), (max(heteros$pcoap$V2) - 0.01), "A", bty="n", cex=1.5)
  text(min(heteros$pcoap$V1), (max(heteros$pcoap$V2) - 0.01), "Bacteria", bty="n", cex=1.5, adj = c(0.2,0))
  box(lwd=2)
  # Add Points & Ellipses
  points(heteros$pcoap$V1, heteros$pcoap$V2,
    pch=19, cex=2.0, bg="gray", col="gray")
  #text(heteros$pcoap$V1, heteros$pcoap$V2,
  #  labels=heteros$pcoap$Treatment, pos=3)
  ordiellipse(cbind(heteros$pcoap$V1, heteros$pcoap$V2),
    heteros$pcoap$Treatment, kind="se", conf=0.95,
    lwd=2, draw = "polygon", col="gray", border = "black", label=TRUE, cex=1)
  #text(-0.14,-0.04,"SRP", cex=1)
  #text(-0.165,0.01, "ATP", cex=1)
  #text(0.01,-0.06, "Mix", cex=1)
  #text(0.11, 0.11, "Ctrl", cex=1)
  #text(0.07, 0.01, "Phyt", cex=1)
  #text(0.11, -0.035, "AEP", cex=1)

# Photos
  x.dim <- c(min(photos$pcoap$V1)-(max(photos$pcoap$V1)*0.1),
    max(photos$pcoap$V1)+(max(photos$pcoap$V1)*0.1))
  y.dim <- c(min(photos$pcoap$V2)-(max(photos$pcoap$V2)*0.1),
    max(photos$pcoap$V2)+(max(photos$pcoap$V2)*0.1))
  # Initiate Plot
  plot(photos$pcoap$V1, photos$pcoap$V2,
    xlab=paste("PCoA 1 (",photos$var1, "%)", sep=""),
    ylab=paste("PCoA 2 (",photos$var2, "%)", sep=""),
    xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",
    xaxt="n", yaxt="n", cex.lab=1.5, cex.axis=1.2, axes=FALSE)
  axis(side = 1, labels=T, lwd.ticks=2, at=c(-0.3,0,0.15), cex.axis=1.2)
  axis(side = 2, labels=T, lwd.ticks=2, at=c(-0.15,0,0.15), cex.axis=1.2)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  #text(min(photos$pcoap$V1), max(photos$pcoap$V2), "B", bty="n", cex=1.5)
  text(min(photos$pcoap$V1), max(photos$pcoap$V2), "Cyanobacteria", bty="n", cex=1.5, adj = c(0.1,0))
  box(lwd=2)
  # Add Points & Ellipses
  points(photos$pcoap$V1, photos$pcoap$V2,
    pch=19, cex=2.0, bg="gray", col="gray")
  #text(photos$pcoap$V1, photos$pcoap$V2,
  #  labels=photos$pcoap$Treatment, pos=3)
  ordiellipse(cbind(photos$pcoap$V1, photos$pcoap$V2),
    photos$pcoap$Treatment, kind="se", conf=0.95,
    lwd=2, draw = "polygon", col="gray", border = "black", label=TRUE, cex=1) 

# Euks
  x.dim <- c(min(euks$pcoap$V1)-(max(euks$pcoap$V1)*0.1),
    max(euks$pcoap$V1)+(max(euks$pcoap$V1)*0.1))
  y.dim <- c(min(euks$pcoap$V2)-(max(euks$pcoap$V2)*0.1),
    max(euks$pcoap$V2)+(max(euks$pcoap$V2)*0.1))
  # Initiate Plot
  plot(euks$pcoap$V1, euks$pcoap$V2,
    xlab=paste("PCoA 1 (",euks$var1, "%)", sep=""),
    ylab=paste("PCoA 2 (",euks$var2, "%)", sep=""),
    xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",
    xaxt="n", yaxt="n", cex.lab=1.5, cex.axis=1.2, axes=FALSE)
  axis(side = 1, labels=T, lwd.ticks=2, at=c(-0.3,0,0.3), cex.axis=1.2)
  axis(side = 2, labels=T, lwd.ticks=2, at=c(-0.3,0,0.3), cex.axis=1.2)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  #text(min(euks$pcoap$V1), max(euks$pcoap$V2), "C", bty="n", cex=1.5)
  text(min(euks$pcoap$V1), max(euks$pcoap$V2), "Algae", bty="n", cex=1.5, adj = c(0.3,0))
  box(lwd=2)
  # Add Points & Ellipses
  points(euks$pcoap$V1, euks$pcoap$V2,
    pch=19, cex=2.0, bg="gray", col="gray")
  #text(euks$pcoap$V1, euks$pcoap$V2,
  #  labels=euks$pcoap$Treatment, pos=3)
  ordiellipse(cbind(euks$pcoap$V1, euks$pcoap$V2),
    euks$pcoap$Treatment, kind="se", conf=0.95,
    lwd=2, draw = "polygon", col="gray", border = "black", label=TRUE, cex=1) 


  # Save Plots
  dev.copy2pdf(file=paste("../graphs/ordinations.pdf",sep=""), useDingbats=FALSE)
  dev.copy(png, file=paste("../graphs/ordinations.png",sep=""), width=72*(4),
    height=72*(10), res=72*2)
  dev.off()