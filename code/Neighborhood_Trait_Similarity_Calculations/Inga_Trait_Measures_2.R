#!/usr/bin/Rscript
library(ape)
library(vegan)
bci_phylo <- read.nexus(file = "./data/BCI_Inga_Phylogeny_PNAS_2017.nex")
plot(bci_phylo)
congen_het_dist <- cophenetic.phylo(bci_phylo)
congen_het_dist <- as.data.frame(congen_het_dist)
congen_het_dist_norm <- congen_het_dist/max(congen_het_dist)
congen_het_sim <- (1-(congen_het_dist/max(congen_het_dist)))


#DEFENSES
#chemsitry

# Sedio's chemistry
chem_sim_sedio  <- read.csv("./data/Trait_Data/species_chemical_similarity_sedio.csv",row.names = 1)
chem_sim_sedio <- (chem_sim_sedio[order(row.names(chem_sim_sedio)), order(names(chem_sim_sedio))])
chem_dist_sedio <- (1-chem_sim_sedio)
chem_dist_sedio_norm <- chem_dist_sedio/max(chem_dist_sedio)
chem_sim_sedio_norm <- (1-chem_dist_sedio_norm)

ants_all <- read.csv("./data/Trait_Data/ants.csv")
ants <- ants_all[c(1:5,7,9:10,13),2:3]
row.names(ants) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
ants_pca <- prcomp(ants,center =T,scale. = T)
ants_dist <-as.data.frame(as.matrix(vegdist(ants_pca$x,method = "euclidean" ,upper = T,diag = T)))
names(ants_dist) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
ants_dist_norm <- ants_dist/max(ants_dist)
ants_sim_norm <- (1-ants_dist_norm)
plot(hclust(as.dist(ants_dist_norm)),main = "Ant Visitation Den")

dev_all <- read.csv("./data/Trait_Data/dev.csv")
dev <- dev_all[c(1:5,7,9:10,13),2:3]
row.names(dev) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
dev_pca <- prcomp(dev,center =T,scale. = T)
dev_dist <-as.data.frame(as.matrix(vegdist(dev_pca$x,method = "euclidean" ,upper = T,diag = T)))
names(dev_dist) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
dev_dist_norm <- (dev_dist/max(dev_dist))
dev_sim_norm <- (1-dev_dist_norm)
plot(hclust(as.dist(dev_dist_norm)),main = "Developement Den")


hair_all <- read.csv("./data/Trait_Data/hair.csv")
hair <- hair_all[c(1:5,7,9:10,13),2:6] 
row.names(hair) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
hair_pca <- prcomp(hair,center =T,scale. = T)
hair_dist <-as.data.frame(as.matrix(vegdist(hair_pca$x,method = "euclidean" ,upper = T,diag = T)))
names(hair_dist) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
hair_dist_norm <- (hair_dist/max(hair_dist))
hair_sim_norm <- (1-hair_dist_norm)
plot(hclust(as.dist(dev_dist_norm)),main = "Developement Den")


#Non_Defenses
wood_all <- read.csv("./data/Trait_Data/wood.csv")
wood  <- wood_all[c(1:5,7,9:10,13),2:3]
row.names(wood) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
wood_pca <- prcomp(wood,center =T,scale. = T)
wood_dist <-as.data.frame(as.matrix(vegdist(wood_pca$x,method = "euclidean" ,upper = T,diag = T)))
names(wood_dist) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
wood_dist_norm <- (wood_dist/max(wood_dist))
wood_sim_norm <- (1-wood_dist_norm)
plot(hclust(as.dist(wood_dist_norm)),main = "Wood_specific_Grav_Dend")

element_all <- read.csv("./data/Trait_Data/element.csv")
element  <- element_all[c(1:5,7,9:10,13),2:13]
row.names(element) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
element_pca <- prcomp(element,center =T,scale. = T)
element_dist <-as.data.frame(as.matrix(vegdist(element_pca$x,method = "euclidean" ,upper = T,diag = T)))
names(element_dist) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
element_dist_norm <- (element_dist/max(element_dist))
element_sim_norm <- (1-element_dist_norm)
plot(hclust(as.dist(element_dist_norm)),main = "Elemental_Composition_Dend")

height_all <- read.csv("./data/Trait_Data/height.csv")
height  <- height_all[c(1:5,7,9:10,13),2:4]
row.names(height) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
height_pca <- prcomp(height,center =T,scale. = T)
height_dist <-as.data.frame(as.matrix(vegdist(height_pca$x,method = "euclidean" ,upper = T,diag = T)))
names(element_dist) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
height_dist_norm <- (height_dist/max(height_dist))
height_sim_norm <- (1-height_dist_norm)
plot(hclust(as.dist(height_dist_norm)),main = "Maximum Size and Height Dend")


leaf_all <- read.csv("./data/Trait_Data/leaf.csv")
leaf  <- leaf_all[c(1:5,7,9:10,13),2:9]
row.names(leaf) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
leaf_pca <- prcomp(leaf,center =T,scale. = T)
leaf_dist <-as.data.frame(as.matrix(vegdist(leaf_pca$x,method = "euclidean" ,upper = T,diag = T)))
names(leaf_dist) <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")
leaf_dist_norm <- (leaf_dist/max(leaf_dist))
leaf_sim_norm <- (1-leaf_dist_norm)
plot(hclust(as.dist(leaf_dist_norm)),main = "Leaf Structure Dend")





