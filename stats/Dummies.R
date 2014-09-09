################################################################################
#                                                                              #
#	Phosphorus Mecocosm Project: Dummy Matrix                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/04/28                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~./GitHub/p-meso/stats")

# Exp Design
design <- read.delim("../data/design.txt")
  design$Treatment <- factor(design$Treatment,
    levels = c('Ctrl', 'SRP', 'AEP', 'ATP', 'Phyt', 'Mix'))

################################################################################
# Calculations and Transformations
################################################################################
  # Dummy Variables #
  treat.no.mix <- design$Treatment
    for(i in 1:length(treat.no.mix)){
      if (treat.no.mix[i] == "Mix") {treat.no.mix[i] = NA}}
  treat.only <- design$Treatment
    for(i in 1:length(treat.only)){
      if (treat.only[i] == "Ctrl") {treat.only[i] = NA}
      else{if (treat.only[i] == "Mix") {treat.only[i] = NA}}}
  ctrl.vs.rest <- design$Treatment == "Ctrl"
    for(i in 1:length(ctrl.vs.rest)){
    if (ctrl.vs.rest[i] == "TRUE") {ctrl.vs.rest[i] = "Ctrl"}
      else (ctrl.vs.rest[i] = "Rest")
    }
  ctrl.vs.rest.wo.mix <- design$Treatment == "Ctrl"
    for(i in 1:length(ctrl.vs.rest.wo.mix)){
      if (ctrl.vs.rest.wo.mix[i] == "TRUE") {ctrl.vs.rest.wo.mix[i] = "Ctrl"}
        else (ctrl.vs.rest.wo.mix[i] = "Rest")}
    for(i in 1:length(design$Treatment)){
      if(design$Treatment[i] == "Mix") {ctrl.vs.rest.wo.mix[i] = NA}
    }
  inor.vs.org <- design$Treatment == "SRP"
    for(i in 1:length(inor.vs.org)){
      if (inor.vs.org[i] == "TRUE") {inor.vs.org[i] = "Inor"}
        else (inor.vs.org[i] = "Org")}
    for(i in 1:length(design$Treatment)){
      if(design$Treatment[i] == "Ctrl"){inor.vs.org[i] = "Ctrl"}
        else{if(design$Treatment[i] == "Mix") inor.vs.org[i] = NA}
    }
  biomolecule <- design$Treatment == "SRP" | design$Treatment == "ATP"
    for(i in 1:length(biomolecule)){
      if (biomolecule[i] == "TRUE") {biomolecule[i] = "Bio"}
        else (biomolecule[i] = "Non-Bio")}
    for(i in 1:length(design$Treatment)){
      if(design$Treatment[i] == "Ctrl"){biomolecule[i] = "Ctrl"}
        else{if(design$Treatment[i] == "Mix") biomolecule[i] = NA}
    }
  complexity.b <- rep(NA, length(design$Treatment))
    for(i in 1:length(design$Treatment)){
      if(design$Treatment[i] == "Mix") {complexity.b[i] = NA}
        else{if(design$Treatment[i] == "Ctrl"){complexity.b[i] = 0}
        else{if(design$Treatment[i] == "SRP"){complexity.b[i]  = 15.2}
        else{if(design$Treatment[i] == "AEP"){complexity.b[i]  = 271.0}
        else{if(design$Treatment[i] == "Phyt"){complexity.b[i] = 457.9}
        else{if(design$Treatment[i] == "ATP"){complexity.b[i]  = 514.1}
    }}}}}}
  complexity.g <- rep(NA, length(design$Treatment))
    for(i in 1:length(design$Treatment)){
      if(design$Treatment[i] == "Mix") {complexity.g[i] = NA}
        else{if(design$Treatment[i] == "Ctrl"){complexity.g[i] = 0}
        else{if(design$Treatment[i] == "SRP"){complexity.g[i]  = -261.97}
        else{if(design$Treatment[i] == "AEP"){complexity.g[i]  = -211.01}
        else{if(design$Treatment[i] == "Phyt"){complexity.g[i] = -1494.96}
        else{if(design$Treatment[i] == "ATP"){complexity.g[i]  = -671.96}
    }}}}}}
  multifunct <- design$Treatment == "Mix"
    for(i in 1:length(multifunct)){
      if (multifunct[i] == "TRUE") {multifunct[i] = "Mix"}
        else (multifunct[i] = "Rest")}
    for(i in 1:length(design$Treatment)){
      if(design$Treatment[i] == "Ctrl") {multifunct[i] = NA}
    }

  # Make Factors
  tnm <- as.factor(treat.no.mix)
  tonly <- as.factor(treat.only)
  cvr <- as.factor(ctrl.vs.rest)
  cvrnm <- as.factor(ctrl.vs.rest.wo.mix)
  io <- as.factor(inor.vs.org)
  bio <- as.factor(biomolecule)
  com.b <- complexity.b
  com.g <- complexity.g
  mul <- as.factor(multifunct)

  # Dummy Matrix
  dummies <- cbind(design, tnm, tonly, cvr, cvrnm, io, bio, com.b, com.g, mul)

  # Write Dummy Matrix To File
  write.table(dummies, file="./p-meso_dummies.txt", sep="\t", row.names=F)