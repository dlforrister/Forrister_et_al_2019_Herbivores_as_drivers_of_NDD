#!/usr/bin/Rscript
library(ape)
library(spatstat)
library(plyr)
library(vegan)
library(lsr)
library(RMySQL)
attach("./data/CTFSRPackage.rdata")
library(ggplot2)
library(reshape2)
library(Hmisc)

# Due the large amount of data generated in this analyses we have employed a SQL database hosted at the University of Utah's Center for High Performance Computing in order to speed up the analysis. Thus, at the end of each iteration in the loop the results are deposited into the aforementioned database. Note, his is not necessary for the analysis, but makes the code much faster. The credentials have been redacted in the published version of this code. As this may limit people's ability to reproduce this analysis we have included CSV exports for all database tables generated in this analysis in ./data/Focal_Tree_Neighborhood/

mydb = dbConnect(MySQL(), user='u6000251', password='********', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')

#fields define the database column names.
fields <- c("`ft_id`, `iteration`, `census`, `sp`, `quadrat`, `gx`, `gy`, `life_event`, `dbh`, `basal_area`, `sc`, `norm_agb`, `norm_grth`, `an`, `tot`, `cs`, `cg`,`cg_def_c_s`, `cg_def_a`, `cg_def_d`, `cg_def_h`, `cg_def_p_pearson',`cg_nondef_w`, `cg_nondef_e`, `cg_nondef_h`, `cg_nondef_l`,`cg_herb_jaccard_leps_saw_ish`,`gap_TL`")

#this reads in info on the focal species that will be used in this analysis. This list determines which species will be analyzed for the remainder of the analysis.
names <- read.csv("./data/Inga_plot_names_complete_traits.csv",head=T)

#Plot data in this analysis are stored as marked point pattern process objects where each tree is stored as a single point in the space defined by the below window which is a rectangle with dimensions 1000 m by 500 m. Each point within this window is then defined by its x and y coordinates in meters from the origin and marked by its plot tag identifier.
owin_plot <<- owin(xrange=c(-1,1001), yrange=c(-1,501))

#This generates 20 equal area anuli surounding a focal tree whose radii vary from 2.81 to 12.6 m.

largest= pi*12.6^2
EEA = largest/20 
ep =EEA/pi
R=12.6
H1 = sqrt(R*R-ep)
an_mid_1 = round((R+H1)/2,1)
R=H1
H2 = sqrt(R*R-ep)
an_mid_2 = round((R+H2)/2,1)
R=H2
H3 = sqrt(R*R-ep)
an_mid_3 = round((R+H3)/2,1)
R=H3
H4 = sqrt(R*R-ep)
an_mid_4 = round((R+H4)/2,1)
R=H4
H5 = sqrt(R*R-ep)
an_mid_5 = round((R+H5)/2,1)
R=H5
H6 = sqrt(R*R-ep)
an_mid_6 = round((R+H6)/2,1)
R=H6
H7 = sqrt(R*R-ep)
an_mid_7 = round((R+H7)/2,1)
R=H7
H8 = sqrt(R*R-ep)
an_mid_8 = round((R+H8)/2,1)
R=H8
H9 = sqrt(R*R-ep)
an_mid_9 = round((R+H9)/2,1)
R=H9
H10 = sqrt(R*R-ep)
an_mid_10 = round((R+H10)/2,1)
R=H10
H11 = sqrt(R*R-ep)
an_mid_11 = round((R+H11)/2,1)
R=H11
H12 = sqrt(R*R-ep)
an_mid_12 = round((R+H12)/2,1)
R=H12
H13 = sqrt(R*R-ep)
an_mid_13 = round((R+H13)/2,1)
R=H13
H14 = sqrt(R*R-ep)
an_mid_14 = round((R+H14)/2,1)
R=H14
H15 = sqrt(R*R-ep)
an_mid_15 = round((R+H15)/2,1)
R=H15
H16 = sqrt(R*R-ep)
an_mid_16 = round((R+H16)/2,1)
R=H16
H17 = sqrt(R*R-ep)
an_mid_17 = round((R+H17)/2,1)
R=H17
H18 = sqrt(R*R-ep)
an_mid_18 = round((R+H18)/2,1)
R=H18
H19 = sqrt(R*R-ep)
an_mid_19 = round((R+H19)/2,1)
an_mid_20 = round((0+H19)/2,1)


Rval <-c(0, 12.6,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,H13,H14,H15,H16,H17,H18,H19)
An <- Rval[order(Rval)]
lab= c(an_mid_1, an_mid_2,  an_mid_3, an_mid_4, an_mid_5,an_mid_6,an_mid_7,  an_mid_8, an_mid_9, an_mid_10, an_mid_11, an_mid_12, an_mid_13, an_mid_14, an_mid_15, an_mid_16, an_mid_17, an_mid_18, an_mid_19,an_mid_20)




#This code reads in Inga trait similarity measurements expluding data on plant phenology.
source("./code/Neighborhood_Trait_Similarity_Calculations/Inga_Trait_Measures_2.R")

#Code for reading in and analyzing data on the phenology of new leaf production
source("./code/Neighborhood_Trait_Similarity_Calculations/Phenology_BCI_in_R.R")

## This reads in the species which we define as gap specialists. See supplementary material for details:
TL_gap <-read.csv("./data/Coley_Kursar_gap_specialists.csv")
TL_gap_sp <- TL_gap$Code


#Code for quantifying herbivore similarity.
source("./code/Neighborhood_Trait_Similarity_Calculations/herbivore_similarity_calculation.R")


#This loads in all neceasary tree census data provided by CTFS. Please see notes in ./data/Tree_Tables in order to access these data from CTFS.

load("./data/Tree_Tables/bci.tree1.rdata")
load("./data/Tree_Tables/bci.tree2.rdata")
load("./data/Tree_Tables/bci.tree3.rdata")
load("./data/Tree_Tables/bci.tree4.rdata")
load("./data/Tree_Tables/bci.tree5.rdata")
load("./data/Tree_Tables/bci.tree6.rdata")
load("./data/Tree_Tables/bci.tree7.rdata")
load("./data/Tree_Tables/bci.tree8.rdata")

bci.full1 <- bci.tree1
bci.full2 <- bci.tree2
bci.full3 <- bci.tree3
bci.full4 <- bci.tree4
bci.full5 <- bci.tree5
bci.full6 <- bci.tree6
bci.full7 <- bci.tree7
bci.full8 <- bci.tree8



