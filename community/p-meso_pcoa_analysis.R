################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment: Community PcoA and PERMANOVA                 #
#   For Heterotrophs, Phototrohops, & Eukaryotes                               #
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
setwd("~/GitHub/p-meso/community/")
source("../functions/p-meso.pcoa.R")

# Heterotrophs
heteros <- pmeso.pcoa(shared = "../mothur/output/tanks.bac.final.shared",
  design = "../data/design.txt", plot.title = "heterotroph_pcoa", level = "0.03")

# Phototrophs
photos <- pmeso.pcoa(shared = "../mothur/output/tanks.photo.final.shared",
  design = "../data/design.txt", plot.title = "phototroph_pcoa", level = "0.03")

# Eukaryotes
euks <- pmeso.pcoa(shared = "../mothur/output/tanks.euk.final.shared",
  design = "../data/design.txt", plot.title = "eukaryote_pcoa", level = "0.03")

