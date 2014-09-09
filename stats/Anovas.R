################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: P-meso Results: Anova                           #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/04/29                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~./GitHub/p-meso/stats")
require("car")
require("reshape")
require("vegan") 
source("../functions/DiversityFunctions.r")
source("../functions/p-meso.diversity.R")

# Import Data 
  # Eco
  eco.data <- read.delim("../data/p-meso_eco.txt", header=T)
  # Chem  
  chem.data <- read.delim("../data/p-meso_chem.txt", header=T)
  # Micro
  bio.data <- read.delim("../data/p-meso_bio.txt", sep="\t", header=T)

  # Exp Design
  design <- read.delim("../data/design.txt")
  design$Treatment <- factor(design$Treatment, 
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))
  # Dummy Matrix
  dummies <- read.delim("p-meso_dummies.txt")
  # Response Matrix
  responses <- merge(merge(chem.data, bio.data, by="Sample"), 
    eco.data, by="Sample")
  # Responses of Interest
  response <- c("DOC.uM", "TN.uM", "TP.uM", "BP.uMC", "BR.uMO2", "BGE", "ChlA2", 
    "Resp_uM", "NPP_uM")  
  # Models of Interest
  model <- colnames(dummies)[2:11]
  
# Calculate Diversity
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
# ANOVAs
################################################################################

# Chemistry ####################################################################
  # Test for a difference in P conc by P-treatment
  p.anova1 <- lm(chem.data$TP.uM ~ dummies$tonly)
  Anova(p.anova1)
  
    p.anova2 <- lm(chem.data$TP.uM ~ dummies$cvr)
  Anova(p.anova2)

  # Test for a difference in N conc by P-treatment
  n.anova1 <- lm(chem.data$TN.uM ~ dummies$tonly)
  Anova(n.anova1)
  
  n.anova2 <- lm(chem.data$TN.uM ~ dummies$cvr)
  Anova(n.anova2)

  # Test for a difference in DOC conc by P-treatment
  doc.anova1 <- lm(chem.data$DOC.uM ~ dummies$tonly)
  Anova(doc.anova1)
  
  doc.anova2 <- lm(chem.data$DOC.uM ~ dummies$cvr)
  Anova(doc.anova2)

  oneway.test(chem.data$DOC.uM ~ dummies$tonly)

# Microbial Physiology #########################################################
  # Test for a difference in BP by P-treatment
  bp.anova1 <- lm(bio.data$BP.uMC ~ dummies$cvr)
  Anova(bp.anova1)

  bp.anova2 <- lm(bio.data$BP.uMC ~ dummies$tonly)
  Anova(bp.anova2)

  # Test for a difference in BR by P-treatment
  br.anova1 <- lm(bio.data$BR.uMO2 ~ dummies$cvr)
  Anova(br.anova1)

  br.anova2 <- lm(bio.data$BR.uMO2 ~ dummies$tonly)
  Anova(br.anova2)

  # Test for a difference in BGE by P-treatment
  bge.anova1 <- lm(bio.data$BGE ~ dummies$cvr)
  Anova(bge.anova1)

  bge.anova2 <- lm(bio.data$BGE ~ dummies$tonly)
  Anova(bge.anova2)

# Food Web and Ecosystem #######################################################
  # ChlA (time 2)
  
  c.anova1 <- oneway.test(eco.data$ChlA2 ~ dummies$cvr)

  
  c.anova1 <- oneway.test(eco.data$ChlA2 ~ dummies$tonly)
  c.anova1

  # Zoop
  z.anovaT1 <- lm(eco.data$Total_Zoop ~ dummies$cvr)
  Anova(z.anovaT1)

  z.anovaT2 <- lm(eco.data$Total_Zoop ~ dummies$tonly)
  Anova(z.anovaT2)

  z.anovaR1 <- lm(eco.data$CladCope ~ dummies$cvr)
  Anova(z.anovaR1)

  z.anovaR2 <- lm(eco.data$CladCope ~ dummies$tonly)
  Anova(z.anovaR2)

  # Respiration
  resp.anova1 <- lm(eco.data$Resp_uM ~ dummies$cvr)
  Anova(resp.anova1)

  resp.anova2 <- lm(eco.data$Resp_uM ~ dummies$tonly)
  Anova(resp.anova2)

  # NPP
  npp.anova1 <- lm(eco.data$NPP_uM ~ dummies$cvr)
  Anova(npp.anova1)

  npp.anova2 <- lm(eco.data$NPP_uM ~ dummies$tonly)
  Anova(npp.anova2)

# Diversity ####################################################################
  # Test for a difference in Bac-Rich by P-treatment
  b.anovaR1 <- lm(heteros$Richness$mean ~ dummies$cvrnm)
  Anova(b.anovaR1)

  b.anovaR2 <- lm(heteros$Richness$mean ~ dummies$tonly)
  Anova(b.anovaR2)

  # Test for a difference in Bac-Even by P-treatment
  b.anovaE1 <- lm(heteros$Evenness$mean ~ dummies$cvrnm)
  Anova(b.anovaE1)

  b.anovaE2 <- lm(heteros$Evenness$mean ~ dummies$tonly)
  Anova(b.anovaE2)

  # Test for a difference in Cyano-Rich by P-treatment
  p.anovaR1 <- lm(photos$Richness$mean ~ dummies$cvrnm)
  Anova(p.anovaR1)

  p.anovaR2 <- lm(photos$Richness$mean ~ dummies$tonly)
  Anova(p.anovaR2)

  # Test for a difference in Cyano-Even by P-treatment
  p.anovaE1 <- lm(photos$Evenness$mean ~ dummies$cvrnm)
  Anova(p.anovaE1)

  p.anovaE2 <- lm(photos$Evenness$mean ~ dummies$tonly)
  Anova(p.anovaE2)

  # Test for a difference in Euk-Rich by P-treatment
  e.anova <- lm(euks$Evenness$mean ~ dummies$Treatment)
  Anova(e.anova)
  
  e.anovaR1 <- lm(euks$Richness$mean ~ dummies$cvrnm)
  Anova(e.anovaR1)
  
  e.anovaR2 <- lm(euks$Richness$mean ~ dummies$tonly)
  Anova(e.anovaR2)

  # Test for a difference in Euk-Even by P-treatment
  e.anovaE1 <- lm(euks$Evenness$mean ~ dummies$cvrnm)
  Anova(e.anovaE1)

  e.anovaE2 <- lm(euks$Evenness$mean ~ dummies$tonly)
  Anova(e.anovaE2)

# END ##########################################################################
