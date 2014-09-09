################################################################################
#                                                                              #
#	Phosphoros Mesocosm Experiment: Microbial Community Indicator Species        #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/27                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/p-meso/community/")
require("reshape")
source("../functions/p-meso.indval.R")

# Species Indicator Values (see Defrene and Legendre 1997)
# Note: Groups - 1 = Ctrl, 2 = SRP, 3 = AEP, 4= ATP, 5= Phyt, 6= Mix
# Heterotrophs
heteros <- pmeso.indval(shared = "../mothur/output/tanks.bac.final.shared",
  design = "../data/design.txt", level = "0.03")

# Phototrophs
photos <- pmeso.indval(shared = "../mothur/output/tanks.photo.final.shared",
  design = "../data/design.txt", level = "0.03")

# Eukaryotes
euks <- pmeso.indval(shared = "../mothur/output/tanks.euk.final.shared",
  design = "../data/design.txt", level = "0.03")

# Combine with Taxonomy
# Heteros
hetero.tax.raw <- (read.delim("../mothur/output/tanks.bac.final.0.03.taxonomy"))
hetero.tax <- transform(hetero.tax.raw, Taxonomy=colsplit(hetero.tax.raw[,3],
 split="\\;", names=c("Domain","Phylum","Class","Order","Family","Genus")))
rownames(hetero.tax) <- hetero.tax[,1]
rownames(heteros) <- gsub("Otu0", "Otu", rownames(heteros)) 
hetero.data <- merge(heteros, hetero.tax, by = "row.names", all.x=T )

# Photos
photo.tax.raw <- (read.delim("../mothur/output/tanks.photo.final.0.03.taxonomy"))
photo.tax <- transform(photo.tax.raw, Taxonomy=colsplit(photo.tax.raw[,3],
 split="\\;", names=c("Domain","Phylum","Class","Order","Family","Genus")))
rownames(photo.tax) <- photo.tax[,1]
photo.data <- merge(photos, photo.tax, by = "row.names", all.x=T )
 
# Euks
euk.tax.raw <- (read.delim("../mothur/output/tanks.euk.final.0.03.taxonomy"))
euk.tax <- transform(euk.tax.raw,Taxonomy=colsplit(euk.tax.raw[,3],
 split="\\;", names=c("Domain","Phylum","Class","Order","Family","Genus")))
rownames(euk.tax) <- euk.tax[,1] 
rownames(euks) <- gsub("Otu0", "Otu", rownames(euks)) 
euk.data <- merge(euks, euk.tax, by = "row.names", all.x=T )


################################################################################

# Heteros grouping SRP and ATP (bio = 1, non-bio =2)
heteros.bio <- pmeso.indval.bio(shared = "../mothur/output/tanks.bac.final.shared", 
 design = "../data/design.txt", dummies = "../stats/p-meso_dummies.txt", 
 level = "0.03")
 
# Combine with Taxonomy
hetero.tax.raw <- (read.delim("../mothur/output/tanks.bac.final.0.03.taxonomy"))
hetero.tax <- transform(hetero.tax.raw, Taxonomy=colsplit(hetero.tax.raw[,3],
 split="\\;", names=c("Domain","Phylum","Class","Order","Family","Genus")))
rownames(hetero.tax) <- hetero.tax[,1]
rownames(heteros.bio) <- gsub("Otu0", "Otu", rownames(heteros.bio)) 
hetero.data <- merge(heteros.bio, hetero.tax, by = "row.names", all.x=T )
hetero.data <- hetero.data[order(hetero.data$group, -hetero.data$indval), ]
hetero.data <- hetero.data[hetero.data$Taxonomy$Phylum != "unclassified(100)",]
