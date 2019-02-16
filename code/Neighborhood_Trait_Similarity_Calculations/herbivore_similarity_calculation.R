#!/usr/bin/Rscript
herb<-read.csv("./data/Trait_Data/inga_herbivore_V3_Feb_27_2017.csv")
library(reshape2)
library(vegan)
library(Hmisc)
library(corrplot)

names(herb) <- c("index","plant.species","plant.sp","collection.number","n","family","genus","species","motu","Order")

#This analysis only includes  lepidoptera and  sawfly herbivores.
herbs <- herb[,c("plant.sp","n","family","species")]
leps_saw <- subset(herbs, herbs$family %nin% c("Chalcidoidea","Curculionidae","Coreidae","Cecidomyiidae","Chrysomelidae","Syrphidae"))


species_leps_saw <- leps_saw[,c(1,2,4)]
species_leps_saw <- species_leps_saw[which(species_leps_saw$plant.sp %in% c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum")),]

species_leps_saw_wide <-dcast(
  species_leps_saw, 
  plant.sp ~ species,
  value.var="n",
  fun.aggregate=sum
)

row.names(species_leps_saw_wide) <- species_leps_saw_wide$plant.sp
species_leps_saw_wide <- species_leps_saw_wide[c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum"),]
species_leps_saw_wide <- species_leps_saw_wide[,-1]
spec_wide_leps_saw_normal <- decostand(species_leps_saw_wide, "total", MARGIN=1)

spe.jaccard_leps_saw <- vegdist(spec_wide_leps_saw_normal, method="jaccard")
herb_dist_jaccard_leps_saw_ish <- as.data.frame(stats:::as.matrix.dist(spe.jaccard_leps_saw))
herb_dist_jaccard_leps_saw_ish_norm <- (herb_dist_jaccard_leps_saw_ish/max(herb_dist_jaccard_leps_saw_ish))
herb_sim_norm_jaccard_leps_saw_ish <- (1- herb_dist_jaccard_leps_saw_ish_norm)




#host dietbreadth figure
species_leps_saw_binary <- species_leps_saw
species_leps_saw_binary$n <- 1

species_leps_saw_binary_wide <-dcast(
  unique(species_leps_saw_binary), 
  plant.sp ~ species,
  value.var="n",
  fun.aggregate=sum
)
species_leps_saw_binary_wide <- species_leps_saw_binary_wide[,-1]
nspp <- data.frame(motu = names(species_leps_saw_binary_wide),nspp = colSums(species_leps_saw_binary_wide))
barplot(table(nspp$nspp))

#Herbivore similarity correlation plot


corrplot::corrplot(as.matrix(herb_sim_norm_jaccard_leps_saw_ish) , is.corr = FALSE , type = "upper",diag = F,method = "color",tl.cex = 3,cl.cex = 1)


