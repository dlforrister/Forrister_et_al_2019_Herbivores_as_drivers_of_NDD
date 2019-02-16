source("./code/Neighborhood_Trait_Similarity_Calculations/Inga_Trait_Measures_2.R")
source("./code/Neighborhood_Trait_Similarity_Calculations/Phenology_BCI_in_R.R")

phylo <- data.frame(t(combn(names(congen_het_sim),2)),dist=t(congen_het_sim)[lower.tri(congen_het_sim)])$dist
ants <- data.frame(t(combn(names(ants_sim_norm),2)),dist=t(ants_sim_norm)[lower.tri(ants_sim_norm)])$dist
chem <- data.frame(t(combn(names(chem_sim_sedio_norm),2)),dist=t(chem_sim_sedio_norm)[lower.tri(chem_sim_sedio_norm)])$dist
dev <- data.frame(t(combn(names(dev_sim_norm),2)),dist=t(dev_sim_norm)[lower.tri(dev_sim_norm)])$dist
hair <- data.frame(t(combn(names(hair_sim_norm),2)),dist=t(hair_sim_norm)[lower.tri(hair_sim_norm)])$dist
phen <- data.frame(t(combn(names(phen_pearson_sim_norm),2)),dist=t(phen_pearson_sim_norm)[lower.tri(phen_pearson_sim_norm)])$dist
wood <- data.frame(t(combn(names(wood_sim_norm),2)),dist=t(wood_sim_norm)[lower.tri(wood_sim_norm)])$dist
element <- data.frame(t(combn(names(element_sim_norm),2)),dist=t(element_sim_norm)[lower.tri(element_sim_norm)])$dist
height <- data.frame(t(combn(names(height_sim_norm),2)),dist=t(height_sim_norm)[lower.tri(height_sim_norm)])$dist 
leaf <- data.frame(t(combn(names(leaf_sim_norm),2)),dist=t(leaf_sim_norm)[lower.tri(leaf_sim_norm)])$dist


triats <- cbind(phylo,ants,chem,dev,hair,phen,wood,element,height,leaf)
library(psych)
pairs.panels(triats,pch=20, cex.cor=1,cex = 3,rug = F,xlim = c(0,1), ylim =c(0,1),ellipses = F, lm = T)
library(ade4)

congen_het_sim <- as.dist(congen_het_sim)
ants_sim_norm <- as.dist(ants_sim_norm)
chem_sim_norm <- as.dist(chem_sim_sedio_norm)
dev_sim_norm <- as.dist(dev_sim_norm)
hair_sim_norm <- as.dist(hair_sim_norm)
phen_pearson_sim_norm <- as.dist(phen_pearson_sim_norm)
wood_sim_norm <- as.dist(wood_sim_norm)
element_sim_norm <- as.dist(element_sim_norm)
height_sim_norm <- as.dist(height_sim_norm)
leaf_sim_norm <- as.dist(leaf_sim_norm)


mantel.rtest(m1 = congen_het_sim, m2 = ants_sim_norm,nrepet = 10000)
mantel.rtest(m1 = congen_het_sim, m2 = chem_sim_norm,nrepet = 10000)
mantel.rtest(m1 = congen_het_sim, m2 = dev_sim_norm,nrepet = 10000)
mantel.rtest(m1 = congen_het_sim, m2 = hair_sim_norm,nrepet = 10000)
mantel.rtest(m1 = congen_het_sim, m2 = phen_pearson_sim_norm,nrepet = 10000)
mantel.rtest(m1 = congen_het_sim, m2 = wood_sim_norm,nrepet = 10000)
mantel.rtest(m1 = congen_het_sim, m2 = element_sim_norm,nrepet = 10000)
mantel.rtest(m1 = congen_het_sim, m2 = height_sim_norm,nrepet = 10000)
mantel.rtest(m1 = congen_het_sim, m2 = leaf_sim_norm,nrepet = 10000)

mantel.rtest(m1 = ants_sim_norm, m2 = chem_sim_norm,nrepet = 10000)
mantel.rtest(m1 = ants_sim_norm, m2 = dev_sim_norm,nrepet = 10000)
mantel.rtest(m1 = ants_sim_norm, m2 = hair_sim_norm,nrepet = 10000)
mantel.rtest(m1 = ants_sim_norm, m2 = phen_pearson_sim_norm,nrepet = 10000)
mantel.rtest(m1 = ants_sim_norm, m2 = wood_sim_norm,nrepet = 10000)
mantel.rtest(m1 = ants_sim_norm, m2 = element_sim_norm,nrepet = 10000)
mantel.rtest(m1 = ants_sim_norm, m2 = height_sim_norm,nrepet = 10000)
mantel.rtest(m1 = ants_sim_norm, m2 = leaf_sim_norm,nrepet = 10000)

mantel.rtest(m1 = chem_sim_norm, m2 = dev_sim_norm,nrepet = 10000)
mantel.rtest(m1 = chem_sim_norm, m2 = hair_sim_norm,nrepet = 10000)
mantel.rtest(m1 = chem_sim_norm, m2 = phen_pearson_sim_norm,nrepet = 10000)
mantel.rtest(m1 = chem_sim_norm, m2 = wood_sim_norm,nrepet = 10000)
mantel.rtest(m1 = chem_sim_norm, m2 = element_sim_norm,nrepet = 10000)
mantel.rtest(m1 = chem_sim_norm, m2 = height_sim_norm,nrepet = 10000)
mantel.rtest(m1 = chem_sim_norm, m2 = leaf_sim_norm,nrepet = 10000)

mantel.rtest(m1 = dev_sim_norm, m2 = hair_sim_norm,nrepet = 10000)
mantel.rtest(m1 = dev_sim_norm, m2 = phen_pearson_sim_norm,nrepet = 10000)
mantel.rtest(m1 = dev_sim_norm, m2 = wood_sim_norm,nrepet = 10000)
mantel.rtest(m1 = dev_sim_norm, m2 = element_sim_norm,nrepet = 10000)
mantel.rtest(m1 = dev_sim_norm, m2 = height_sim_norm,nrepet = 10000)
mantel.rtest(m1 = dev_sim_norm, m2 = leaf_sim_norm,nrepet = 10000)

mantel.rtest(m1 = hair_sim_norm, m2 = phen_pearson_sim_norm,nrepet = 10000)
mantel.rtest(m1 = hair_sim_norm, m2 = wood_sim_norm,nrepet = 10000)
mantel.rtest(m1 = hair_sim_norm, m2 = element_sim_norm,nrepet = 10000)
mantel.rtest(m1 = hair_sim_norm, m2 = height_sim_norm,nrepet = 10000)
mantel.rtest(m1 = hair_sim_norm, m2 = leaf_sim_norm,nrepet = 10000)

mantel.rtest(m1 = phen_pearson_sim_norm, m2 = wood_sim_norm,nrepet = 10000)
mantel.rtest(m1 = phen_pearson_sim_norm, m2 = element_sim_norm,nrepet = 10000)
mantel.rtest(m1 = phen_pearson_sim_norm, m2 = height_sim_norm,nrepet = 10000)
mantel.rtest(m1 = phen_pearson_sim_norm, m2 = leaf_sim_norm,nrepet = 10000)

mantel.rtest(m1 = wood_sim_norm, m2 = element_sim_norm,nrepet = 10000)
mantel.rtest(m1 = wood_sim_norm, m2 = height_sim_norm,nrepet = 10000)
mantel.rtest(m1 = wood_sim_norm, m2 = leaf_sim_norm,nrepet = 10000)

mantel.rtest(m1 = element_sim_norm, m2 = height_sim_norm,nrepet = 10000)
mantel.rtest(m1 = element_sim_norm, m2 = leaf_sim_norm,nrepet = 10000)

mantel.rtest(m1 = height_sim_norm, m2 = leaf_sim_norm,nrepet = 10000)
