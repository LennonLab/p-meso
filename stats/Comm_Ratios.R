################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment:  Microbial Community Ratios                  #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/26                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~./GitHub/p-meso/stats/")
source("../functions/DiversityFunctions.r")
require("vegan")
se <- function(x, ...){sd(x, ...)/sqrt(length(x))}

# Exp Design
design <- read.delim("../data/design.txt")
  design$Treatment <- factor(design$Treatment,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

# Calculate the number of observations for each functional group 
heteros <- count.groups(read.otu(shared = "../mothur/output/tanks.bac.final.shared", 
  cutoff = "0.03"))

cyanos <- count.groups(read.otu(shared = "../mothur/output/tanks.photo.final.shared",
  cutoff = "0.03"))
  
ratios <- cbind(design, heteros, cyanos)[-1]
ratios$CyanoHetero <- ratios$cyanos / ratios$heteros


ratio.means <- tapply(ratios$CyanoHetero, ratios$Treatment, mean)
ratio.se <- tapply(ratios$CyanoHetero, ratios$Treatment, se)


  # Dummy Matrix
  dummies <- read.delim("p-meso_dummies.txt")
ratio.lm <- lm(ratios$CyanoHetero ~ dummies$cvrnm)
Anova(ratio.lm)
  



