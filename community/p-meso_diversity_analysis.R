################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment: Diversity                                    #
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
source("../functions/p-meso.diversity.R")

# Heterotrophs
heteros <- pmeso.diversity(shared = "../mothur/output/tanks.bac.final.shared",
  design = "../data/design.txt", plot.title = "heterotroph_div", level = "0.03",
  size = 10000)

# Phototrophs
photos <- pmeso.diversity(shared = "../mothur/output/tanks.photo.final.shared",
  design = "../data/design.txt", plot.title = "phototroph_div", level = "0.03",
  size = 10000)

# Eukaryotes
euks <- pmeso.diversity(shared = "../mothur/output/tanks.euk.final.shared",
  design = "../data/design.txt", plot.title = "eukaryote_div", level = "unique",
  size = 100)

