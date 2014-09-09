################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment:  Microbial Community Diversity               #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/26                                                      #
#                                                                              #
################################################################################

pmeso.diversity <- function(shared = " ", design = " ", plot.title = "out",
  level = "0.03", size = 10000){

  # Setup Work Environment
  source("../functions/DiversityFunctions.r")
  require("vegan")
  se <- function(x, ...){sd(x, ...)/sqrt(length(x))}
  design <- read.delim(design, header=T, row.names=1)
  
# Alpha Diversity with Resampling
  rich <- round(richness.iter(input = shared, cutoff = level, size = size,
    iters = 100), 3)
  even <- round(evenness.iter(input = shared, cutoff = level, size = size,
    iters = 100), 3)
  shan <- round(diversity.iter(input = shared, index = "shannon",
    cutoff = level, size = size, iters = 100), 3)

  # Richness Summary
  rich_data <- merge(design, rich, by="row.names")
  rich_data$mean <- round(apply(rich_data[3:26], 1, mean, na.rm = TRUE),3)
  rich_data$se <- round(apply(rich_data[3:26], 1, se, na.rm = TRUE), 3)
  rich_data$Design <- design$Treatment
  rich_data$Design <- factor(rich_data$Design,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Evenness Summary
  even_data <- merge(design, even, by="row.names")
  even_data$mean <- round(apply(even_data[3:26], 1, mean, na.rm = TRUE),3)
  even_data$se <- round(apply(even_data[3:26], 1, se, na.rm = TRUE),3)
  even_data$Design <- design$Treatment
  even_data$Design <- factor(even_data$Design,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Shannon Diversity Summary
  shan_data <- merge(design, shan, by="row.names")
  shan_data$mean <- round(apply(shan_data[3:26], 1, mean, na.rm = TRUE),3)
  shan_data$se <- round(apply(shan_data[3:26], 1, se, na.rm = TRUE),3)
  shan_data$Design <- design$Treatment
  shan_data$Design <- factor(shan_data$Design,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Plot Parameters
  windows.options(width=6, height=12)
  par(mfrow=c(3,1), mar=c(0.25,5,0.25,1), oma=c(3,1,1,1)+0.1 )

  # Richness Plot A
  rich.means <- tapply(rich_data$mean, rich_data$Design, mean)
  rich.se <- tapply(rich_data$mean, rich_data$Design, se)
  dim.r <- round(c(min(rich.means)-(max(rich.se)*1.1),
    max(rich.means)+(max(rich.se)*1.1)),2)
  plot(rich.means, ylim=dim.r, xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  arrows(x0 = 1:6, y0 = rich.means, y1 = rich.means - rich.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = rich.means, y1 = rich.means + rich.se, angle = 90,
    length=0.05, lwd = 2)
  title(ylab = "Taxa Richness", cex.lab = 2, line = 3.5)
  box(lwd=2)
  axis(side = 1, labels=F, lwd.ticks=2)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

  # Evenness Plot A
  even.means <- tapply(even_data$mean, even_data$Design, mean)
  even.se <- tapply(even_data$mean, even_data$Design, se)
  dim.e <- round(c(min(even.means)-(max(even.se)*1.1),
    max(even.means)+(max(even.se)*1.1)),2)
  plot(even.means, ylim=dim.e, xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  arrows(x0 = 1:6, y0 = even.means, y1 = even.means - even.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = even.means, y1 = even.means + even.se, angle = 90,
    length=0.05, lwd = 2)
  title(ylab = "Taxa Evenness", cex.lab = 2, line = 3.5)
  box(lwd=2)
  axis(side = 1, labels=F, lwd.ticks=2)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

  # Shannon Plot A
  shan.means <- tapply(shan_data$mean, shan_data$Design, mean)
  shan.se <- tapply(shan_data$mean, shan_data$Design, se)
  dim.s <- round(c(min(shan.means)-(max(shan.se)*1.1),
    max(shan.means)+(max(shan.se)*1.1)),2)
  plot(shan.means, ylim=dim.s, xlim=c(0.5,6.5), pch=15, cex=2, yaxt="n",
    cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
  points(shan.means, pch=15, cex=2)
  arrows(x0 = 1:6, y0 = shan.means, y1 = shan.means - shan.se, angle = 90,
    length=0.05, lwd = 2)
  arrows(x0 = 1:6, y0 = shan.means, y1 = shan.means + shan.se, angle = 90,
    length=0.05, lwd = 2)
  title(ylab = "Shannon Diversity", cex.lab = 2, line = 3.5)
  box(lwd=2)
  axis(side = 1, labels=c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'), at=1:6,
    lwd.ticks=2, cex.axis=1.25)
  axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
  axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
  axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

  # Save Plots
  dev.copy2pdf(file=paste("../graphs/",plot.title,".pdf",sep=""))
  dev.copy(png, file=paste("../graphs/",plot.title,".png",sep=""),
    bg="white", width=72*(4*4), height=72*(8*4), res=72*4)
  dev.off()
  windows.options(reset=TRUE)
  
  counts <- count.groups(read.otu(shared, cutoff = level))
  coverage <- coverage(input= shared, cutoff = level, size = size)
  output <- list("total.seq.coverage" = counts, 
    "subsampling.goods.coverage" = coverage)
  return(output)
  }
  
pmeso.diversity.calc <- function(shared = " ", design = " ", level = "0.03", 
  size = 10000){

  # Setup Work Environment
  source("../functions/DiversityFunctions.r")
  require("vegan")
  se <- function(x, ...){sd(x, ...)/sqrt(length(x))}
  design <- read.delim(design, header=T, row.names=1)
  
# Alpha Diversity with Resampling
  rich <- round(richness.iter(input = shared, cutoff = level, size = size,
    iters = 100), 3)
  even <- round(evenness.iter(input = shared, cutoff = level, size = size,
    iters = 100), 3)
  shan <- round(diversity.iter(input = shared, index = "shannon",
    cutoff = level, size = size, iters = 100), 3)

  # Richness Summary
  rich_data <- merge(design, rich, by="row.names")
  rich_data$mean <- round(apply(rich_data[3:26], 1, mean, na.rm = TRUE),3)
  rich_data$se <- round(apply(rich_data[3:26], 1, se, na.rm = TRUE), 3)
  rich_data$Design <- design$Treatment
  rich_data$Design <- factor(rich_data$Design,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Evenness Summary
  even_data <- merge(design, even, by="row.names")
  even_data$mean <- round(apply(even_data[3:26], 1, mean, na.rm = TRUE),3)
  even_data$se <- round(apply(even_data[3:26], 1, se, na.rm = TRUE),3)
  even_data$Design <- design$Treatment
  even_data$Design <- factor(even_data$Design,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Shannon Diversity Summary
  shan_data <- merge(design, shan, by="row.names")
  shan_data$mean <- round(apply(shan_data[3:26], 1, mean, na.rm = TRUE),3)
  shan_data$se <- round(apply(shan_data[3:26], 1, se, na.rm = TRUE),3)
  shan_data$Design <- design$Treatment
  shan_data$Design <- factor(shan_data$Design,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

  # Richness
  rich.means <- tapply(rich_data$mean, rich_data$Design, mean)
  rich.se <- tapply(rich_data$mean, rich_data$Design, se)

  # Evenness
  even.means <- tapply(even_data$mean, even_data$Design, mean)
  even.se <- tapply(even_data$mean, even_data$Design, se)

  # Shannon
  shan.means <- tapply(shan_data$mean, shan_data$Design, mean)
  shan.se <- tapply(shan_data$mean, shan_data$Design, se)
  
  return(list("Richness" = rich_data, "Evenness" = even_data, 
    "Shannon" = shan_data))
  }