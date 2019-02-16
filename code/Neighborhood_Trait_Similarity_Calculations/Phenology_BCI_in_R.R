library(circular)
library(lubridate)
library(tidyr)
library(vegan)

phen_raw <- read.csv("./data/Trait_Data/Phenology_Raw.csv")

phen_pearson_comp <- data.frame(matrix(nrow = length(unique(phen_raw$SPP)),ncol = length(unique(phen_raw$SPP))))
row.names(phen_pearson_comp) <- levels(phen_raw$SPP)
names(phen_pearson_comp) <- levels(phen_raw$SPP)

phen_simple <- phen_raw[,c("SPP","Angle","PERCENT")]
sp.comb <- combn(unique(phen_simple$SPP), 2, simplify=T)

for (f in 1:ncol(sp.comb)) {
single_sp_phen_1 <- subset(phen_simple, SPP == sp.comb[1,f])
temp_norm_1 <- single_sp_phen_1$PERCENT[1:10]/sum(single_sp_phen_1$PERCENT[1:10])*100
floor_norm_1 <- floor(temp_norm_1)
floor_norm_1[order(temp_norm_1-floor_norm_1,decreasing = T)[1:(100 - sum(floor_norm_1))]] <- floor_norm_1[order(temp_norm_1-floor_norm_1,decreasing = T)[1:(100 - sum(floor_norm_1))]]+1
single_sp_phen_1$PERCENT[1:10] <- floor_norm_1


temp_norm_2 <- single_sp_phen_1[11:22,]$PERCENT <- single_sp_phen_1$PERCENT[11:22]/sum(single_sp_phen_1$PERCENT[11:22])*100
floor_norm_2 <- floor(temp_norm_2)
floor_norm_2[order(temp_norm_2-floor_norm_2,decreasing = T)[1:(100 - sum(floor_norm_2))]] <- floor_norm_2[order(temp_norm_2-floor_norm_2,decreasing = T)[1:(100 - sum(floor_norm_2))]]+1
single_sp_phen_1$PERCENT[11:22] <- floor_norm_2

temp_norm_3 <- single_sp_phen_1[23:34,]$PERCENT <- single_sp_phen_1$PERCENT[23:34]/sum(single_sp_phen_1$PERCENT[23:34])*100
floor_norm_3 <- floor(temp_norm_2)
floor_norm_3[order(temp_norm_3-floor_norm_3,decreasing = T)[1:(100 - sum(floor_norm_3))]] <- floor_norm_3[order(temp_norm_3-floor_norm_3,decreasing = T)[1:(100 - sum(floor_norm_3))]]+1
single_sp_phen_1$PERCENT[23:34] <- floor_norm_3

temp_norm_4 <-single_sp_phen_1[35:45,]$PERCENT <- single_sp_phen_1$PERCENT[35:45]/sum(single_sp_phen_1$PERCENT[35:45])*100
floor_norm_4 <- floor(temp_norm_4)
floor_norm_4[order(temp_norm_4-floor_norm_4,decreasing = T)[1:(100 - sum(floor_norm_4))]] <- floor_norm_4[order(temp_norm_4-floor_norm_4,decreasing = T)[1:(100 - sum(floor_norm_4))]]+1
single_sp_phen_1$PERCENT[35:45] <- floor_norm_4

angles_1 <- rep(single_sp_phen_1[,2], single_sp_phen_1[,3])
length(angles_1)
angles.circular_1 <- circular(angles_1,type = "angles",units = "degrees",template = "geographic", zero = 0,modulo= "2pi", rotation = "counter")

single_sp_phen_2 <- subset(phen_simple, SPP == sp.comb[2,f])
temp_norm_1 <- single_sp_phen_2$PERCENT[1:10]/sum(single_sp_phen_2$PERCENT[1:10])*100
floor_norm_1 <- floor(temp_norm_1)
floor_norm_1[order(temp_norm_1-floor_norm_1,decreasing = T)[1:(100 - sum(floor_norm_1))]] <- floor_norm_1[order(temp_norm_1-floor_norm_1,decreasing = T)[1:(100 - sum(floor_norm_1))]]+1
single_sp_phen_2$PERCENT[1:10] <- floor_norm_1


temp_norm_2 <- single_sp_phen_2[11:22,]$PERCENT <- single_sp_phen_2$PERCENT[11:22]/sum(single_sp_phen_2$PERCENT[11:22])*100
floor_norm_2 <- floor(temp_norm_2)
floor_norm_2[order(temp_norm_2-floor_norm_2,decreasing = T)[1:(100 - sum(floor_norm_2))]] <- floor_norm_2[order(temp_norm_2-floor_norm_2,decreasing = T)[1:(100 - sum(floor_norm_2))]]+1
single_sp_phen_2$PERCENT[11:22] <- floor_norm_2

temp_norm_3 <- single_sp_phen_2[23:34,]$PERCENT <- single_sp_phen_2$PERCENT[23:34]/sum(single_sp_phen_2$PERCENT[23:34])*100
floor_norm_3 <- floor(temp_norm_2)
floor_norm_3[order(temp_norm_3-floor_norm_3,decreasing = T)[1:(100 - sum(floor_norm_3))]] <- floor_norm_3[order(temp_norm_3-floor_norm_3,decreasing = T)[1:(100 - sum(floor_norm_3))]]+1
single_sp_phen_2$PERCENT[23:34] <- floor_norm_3

temp_norm_4 <-single_sp_phen_2[35:45,]$PERCENT <- single_sp_phen_2$PERCENT[35:45]/sum(single_sp_phen_2$PERCENT[35:45])*100
floor_norm_4 <- floor(temp_norm_4)
floor_norm_4[order(temp_norm_4-floor_norm_4,decreasing = T)[1:(100 - sum(floor_norm_4))]] <- floor_norm_4[order(temp_norm_4-floor_norm_4,decreasing = T)[1:(100 - sum(floor_norm_4))]]+1
single_sp_phen_2$PERCENT[35:45] <- floor_norm_4

angles_2 <- rep(single_sp_phen_2[,2], single_sp_phen_2[,3])
length(angles_2)
angles.circular_2 <- circular(angles_2,type = "angles",units = "degrees",template = "geographic", zero = 0,modulo= "2pi", rotation = "counter")

pearson_cor <- cor.circular(angles.circular_1,angles.circular_2)
phen_pearson_comp[as.character(sp.comb[2,f]),as.character(sp.comb[1,f])] <-  pearson_cor
}

for(i in 1:nrow(phen_pearson_comp)) {
  for(j in i:nrow(phen_pearson_comp)) {
    if(i==j) phen_pearson_comp[i,j] <- 1 else 
    if(is.na(phen_pearson_comp[j,i])) phen_pearson_comp[j,i] <- phen_pearson_comp[i,j] else 
      phen_pearson_comp[i,j] <- phen_pearson_comp[j,i]
  }
}
phen_pearson_comp <- abs(phen_pearson_comp)
phen_pearson_dist <- 1- phen_pearson_comp
phen_pearson_dist_norm <- phen_pearson_dist/max(phen_pearson_dist)
phen_pearson_sim_norm <- (1-phen_pearson_dist_norm)