#Code for normalizing focal tree similarity
#Focal trees Neighbohood files were generated using the Plot_Anal_2018_Nov_13 Rscript which was run on the Center for high performance computing cluster at the University of Utah and deposited into a sequal database in order to decrease computation time. We have downloaded this SQL database from the plot analysis as a csv file in order to simplify the rest of the analysis and figure generation.


library(RMySQL)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggplot2)

Null_Growth_Survival <- function(fixed_null){
Null_Fixed <- read.csv("./data/Focal_Tree_Neighborhood/10_Meter_Analysis/BCI_INGA_Null_fixed_2018_Nov_13.csv",row.names = 1)
Null <- read.csv("./data/Focal_Tree_Neighborhood/10_Meter_Analysis/BCI_INGA_Null_2018_Nov_13.csv",row.names = 1)
BCI <- read.csv("./data/Focal_Tree_Neighborhood/10_Meter_Analysis/BCI_INGA_2018_Nov_13.csv")
BCI <- BCI[,-1]

if(fixed_null == T) {Null <- Null_Fixed}

Null <- subset(Null,an <= 10) # check this parameter later
Null_dist <- Null[,c(1:13)]
Null_dist$tot <- Null$tot/Null$an
Null_dist$cs <- Null$cs/Null$an
Null_dist$cg <- Null$cg/Null$an


Null_dist$cg_def_c_s <- Null$cg_def_c_s/Null$an
Null_dist$cg_def_a <- Null$cg_def_a/Null$an
Null_dist$cg_def_d <- Null$cg_def_d/Null$an
Null_dist$cg_def_h <- Null$cg_def_h/Null$an
Null_dist$cg_def_p_pearson <- Null$cg_def_p_pearson/Null$an

Null_dist$cg_nondef_w <- Null$cg_nondef_w/Null$an
Null_dist$cg_nondef_e <- Null$cg_nondef_e/Null$an
Null_dist$cg_nondef_h <- Null$cg_nondef_h/Null$an
Null_dist$cg_nondef_l <- Null$cg_nondef_l/Null$an


Null_dist$cg_herb_jaccard_leps_saw_ish <- Null$cg_herb_jaccard_leps_saw_ish/Null$an

Null_dist$gap_TL <- Null$gap_TL
Null_dist$id_cen <- paste(Null_dist$ft_id,Null_dist$census,sep="_")

Null_dist_2 <- Null_dist[,c(ncol(Null_dist),4:13)]
Null_dist_2 <- Null_dist_2[order(Null_dist_2$id_cen),]
Null_dist_2 <- unique(Null_dist_2)


Null_dist_2$tot <- aggregate((Null_dist$tot), by=list(Null_dist$id_cen), FUN=sum)$x
Null_dist_2$cs <- aggregate((Null_dist$cs), by=list(Null_dist$id_cen), FUN=sum)$x
Null_dist_2$cg <- aggregate((Null_dist$cg), by=list(Null_dist$id_cen), FUN=sum)$x

Null_dist_2$cg_def_c_s <- aggregate((Null_dist$cg_def_c_s), by=list(Null_dist$id_cen), FUN=sum)$x#with DBH weighting
Null_dist_2$cg_def_a <- aggregate((Null_dist$cg_def_a), by=list(Null_dist$id_cen), FUN=sum)$x
Null_dist_2$cg_def_d <- aggregate((Null_dist$cg_def_d), by=list(Null_dist$id_cen), FUN=sum)$x
Null_dist_2$cg_def_h <- aggregate((Null_dist$cg_def_h), by=list(Null_dist$id_cen), FUN=sum)$x

Null_dist_2$cg_def_p_pearson <- aggregate((Null_dist$cg_def_p_pearson), by=list(Null_dist$id_cen), FUN=sum)$x

Null_dist_2$cg_nondef_w <- aggregate((Null_dist$cg_nondef_w), by=list(Null_dist$id_cen), FUN=sum)$x
Null_dist_2$cg_nondef_e <- aggregate((Null_dist$cg_nondef_e), by=list(Null_dist$id_cen), FUN=sum)$x
Null_dist_2$cg_nondef_h <- aggregate((Null_dist$cg_nondef_h), by=list(Null_dist$id_cen), FUN=sum)$x
Null_dist_2$cg_nondef_l <- aggregate((Null_dist$cg_nondef_l), by=list(Null_dist$id_cen), FUN=sum)$x


Null_dist_2$cg_herb_jaccard_leps_saw_ish <- aggregate((Null_dist$cg_herb_jaccard_leps_saw_ish), by=list(Null_dist$id_cen), FUN=sum)$x

Null_dist_2$gap_TL <- aggregate((Null_dist$gap_TL), by=list(Null_dist$id_cen), FUN=mean)$x

Null_sum <-data.frame(tot_mean = mean(Null_dist_2$tot))
Null_sum$tot_sd <- sd(Null_dist_2$tot)
Null_sum$cs_mean <- mean(Null_dist_2$cs)
Null_sum$cs_sd <- sd(Null_dist_2$cs)
Null_sum$cg_mean <- mean(Null_dist_2$cg)
Null_sum$cg_sd <- sd(Null_dist_2$cg)

Null_sum$cg_def_c_s_mean <- mean(Null_dist_2$cg_def_c_s)
Null_sum$cg_def_c_s_sd <- sd(Null_dist_2$cg_def_c_s)


Null_sum$cg_def_a_mean <- mean(Null_dist_2$cg_def_a)
Null_sum$cg_def_a_sd <- sd(Null_dist_2$cg_def_a)
Null_sum$cg_def_d_mean <- mean(Null_dist_2$cg_def_d)
Null_sum$cg_def_d_sd <- sd(Null_dist_2$cg_def_d)

Null_sum$cg_def_h_mean <- mean(Null_dist_2$cg_def_h)
Null_sum$cg_def_h_sd <- sd(Null_dist_2$cg_def_h)

Null_sum$cg_def_p_pearson_mean <- mean(Null_dist_2$cg_def_p_pearson)
Null_sum$cg_def_p_pearson_sd <- sd(Null_dist_2$cg_def_p_pearson)

Null_sum$cg_nondef_w_mean <- mean(Null_dist_2$cg_nondef_w)
Null_sum$cg_nondef_w_sd <- sd(Null_dist_2$cg_nondef_w)

Null_sum$cg_nondef_e_mean <- mean(Null_dist_2$cg_nondef_e)
Null_sum$cg_nondef_e_sd <- sd(Null_dist_2$cg_nondef_e)

Null_sum$cg_nondef_h_mean <- mean(Null_dist_2$cg_nondef_h)
Null_sum$cg_nondef_h_sd <- sd(Null_dist_2$cg_nondef_h)

Null_sum$cg_nondef_l_mean <- mean(Null_dist_2$cg_nondef_l)
Null_sum$cg_nondef_l_sd <- sd(Null_dist_2$cg_nondef_l)



Null_sum$cg_herb_jaccard_leps_saw_ish_mean <- mean(Null_dist_2$cg_herb_jaccard_leps_saw_ish)
Null_sum$cg_herb_jaccard_leps_saw_ish_sd <- sd(Null_dist_2$cg_herb_jaccard_leps_saw_ish)


Null_sum$gap_TL_mean <- mean(Null_dist_2$gap_TL)
Null_sum$gap_TL_sd <- sd(Null_dist_2$gap_TL)


BCI <- subset(BCI,an <= 10)
BCI_dist <- BCI[,c(1:13)]

BCI_dist$tot <- BCI$tot/BCI$an
BCI_dist$cs <- BCI$cs/BCI$an
BCI_dist$cg <- BCI$cg/BCI$an
BCI_dist$cg_def_c_s <- BCI$cg_def_c_s/BCI$an
BCI_dist$cg_def_a <- BCI$cg_def_a/BCI$an
BCI_dist$cg_def_d <- BCI$cg_def_d/BCI$an
BCI_dist$cg_def_h <- BCI$cg_def_h/BCI$an
BCI_dist$cg_def_p_pearson <- BCI$cg_def_p_pearson/BCI$an

BCI_dist$cg_nondef_w <- BCI$cg_nondef_w/BCI$an
BCI_dist$cg_nondef_e <- BCI$cg_nondef_e/BCI$an
BCI_dist$cg_nondef_h <- BCI$cg_nondef_h/BCI$an
BCI_dist$cg_nondef_l <- BCI$cg_nondef_l/BCI$an


BCI_dist$cg_herb_jaccard_leps_saw_ish <- BCI$cg_herb_jaccard_leps_saw_ish/BCI$an
BCI_dist$gap_TL <- BCI$gap_TL
BCI_dist$id_cen <- paste(BCI_dist$ft_id,BCI_dist$census,sep="_")

BCI_dist_2 <- BCI_dist[,c(ncol(BCI_dist),4:13)]
BCI_dist_2 <- BCI_dist_2[order(BCI_dist_2$id_cen),]
BCI_dist_2 <- unique(BCI_dist_2)

BCI_dist_2$tot <- aggregate((BCI_dist$tot), by=list(BCI_dist$id_cen), FUN=sum)$x
BCI_dist_2$cs <- aggregate((BCI_dist$cs), by=list(BCI_dist$id_cen), FUN=sum)$x
BCI_dist_2$cg <- aggregate((BCI_dist$cg), by=list(BCI_dist$id_cen), FUN=sum)$x

BCI_dist_2$cg_def_c_s <- aggregate((BCI_dist$cg_def_c_s), by=list(BCI_dist$id_cen), FUN=sum)$x#with DBH weighting
BCI_dist_2$cg_def_a <- aggregate((BCI_dist$cg_def_a), by=list(BCI_dist$id_cen), FUN=sum)$x
BCI_dist_2$cg_def_d <- aggregate((BCI_dist$cg_def_d), by=list(BCI_dist$id_cen), FUN=sum)$x
BCI_dist_2$cg_def_h <- aggregate((BCI_dist$cg_def_h), by=list(BCI_dist$id_cen), FUN=sum)$x
BCI_dist_2$cg_def_p_pearson <- aggregate((BCI_dist$cg_def_p_pearson), by=list(BCI_dist$id_cen), FUN=sum)$x

BCI_dist_2$cg_nondef_w <- aggregate((BCI_dist$cg_nondef_w), by=list(BCI_dist$id_cen), FUN=sum)$x
BCI_dist_2$cg_nondef_e <- aggregate((BCI_dist$cg_nondef_e), by=list(BCI_dist$id_cen), FUN=sum)$x
BCI_dist_2$cg_nondef_h <- aggregate((BCI_dist$cg_nondef_h), by=list(BCI_dist$id_cen), FUN=sum)$x
BCI_dist_2$cg_nondef_l <- aggregate((BCI_dist$cg_nondef_l), by=list(BCI_dist$id_cen), FUN=sum)$x


BCI_dist_2$cg_herb_jaccard_leps_saw_ish <- aggregate((BCI_dist$cg_herb_jaccard_leps_saw_ish), by=list(BCI_dist$id_cen), FUN=sum)$x

BCI_dist_2$gap_TL <- aggregate((BCI_dist$gap_TL), by=list(BCI_dist$id_cen), FUN=mean)$x


BCI_Norm <- BCI_dist_2[,c(1:11)]

BCI_Norm$tot <- (BCI_dist_2$tot - Null_sum$tot_mean)/Null_sum$tot_sd
BCI_Norm$cs <- (BCI_dist_2$cs - Null_sum$cs_mean)/Null_sum$cs_sd
BCI_Norm$cg <- (BCI_dist_2$cg - Null_sum$cg_mean)/Null_sum$cg_sd

BCI_Norm$cg_def_c_s <- (BCI_dist_2$cg_def_c_s - Null_sum$cg_def_c_s_mean)/Null_sum$cg_def_c_s_sd

BCI_Norm$cg_def_a <- (BCI_dist_2$cg_def_a - Null_sum$cg_def_a_mean)/Null_sum$cg_def_a_sd
BCI_Norm$cg_def_d <- (BCI_dist_2$cg_def_d - Null_sum$cg_def_d_mean)/Null_sum$cg_def_d_sd
BCI_Norm$cg_def_h <- (BCI_dist_2$cg_def_h - Null_sum$cg_def_h_mean)/Null_sum$cg_def_h_sd
BCI_Norm$cg_def_p_pearson <- (BCI_dist_2$cg_def_p_pearson - Null_sum$cg_def_p_pearson_mean)/Null_sum$cg_def_p_pearson_sd

BCI_Norm$cg_nondef_w <- (BCI_dist_2$cg_nondef_w - Null_sum$cg_nondef_w_mean)/Null_sum$cg_nondef_w_sd
BCI_Norm$cg_nondef_e <- (BCI_dist_2$cg_nondef_e - Null_sum$cg_nondef_e_mean)/Null_sum$cg_nondef_e_sd
BCI_Norm$cg_nondef_h <- (BCI_dist_2$cg_nondef_h - Null_sum$cg_nondef_h_mean)/Null_sum$cg_nondef_h_sd
BCI_Norm$cg_nondef_l <- (BCI_dist_2$cg_nondef_l - Null_sum$cg_nondef_l_mean)/Null_sum$cg_nondef_l_sd


BCI_Norm$cg_herb_jaccard_leps_saw_ish <- (BCI_dist_2$cg_herb_jaccard_leps_saw_ish - Null_sum$cg_herb_jaccard_leps_saw_ish_mean)/Null_sum$cg_herb_jaccard_leps_saw_ish_sd

BCI_Norm$gap_TL <- (BCI_dist_2$gap_TL - Null_sum$gap_TL_mean)/Null_sum$gap_TL_sd
BCI_Norm <<- BCI_Norm
}



