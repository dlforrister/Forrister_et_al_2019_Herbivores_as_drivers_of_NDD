# About
This repository contains data as well as the R code used in the data analysis of the 2019 Forrister et al. paper entitled 'Herbivores as drivers of negative density dependence in tropical forest saplings'.

Dependencies:

R code operates using R version 3.4.3 and requires packages, foreach, doParallel, RMySQL, circular, lubridate, vegan, ape, reshape2, Hmisc, corrplot, spatstat, plyr, lsr, ggplot2, picante, tidyr, lme4, lmerTest, MuMIn, splitstackshape, rgdal, spdep, psych, ade4.

Additionally, the CTFS R Package is required. This can be acquired at http://ctfs.si.edu/Public/CTFSRPackage/.

Structure:

All data used in the analysis are available in the data folder, with the exception of the Barro Colorado Forest Census Plot Data which is archived and openly available for download online from the CTFS & ForestGeo http://ctfs.si.edu/datarequest/. 

The code is broken into two sections. The first (Neighborhood_Trait_Similarity_Calculations) deals with calculating neighborhood similarity for all focal trees. This section of the code relies heavily on a SQL database that is hosted by the Center For High Performance Computing at the University of Utah. To allow easier data access all data generated during this first stage is available as csv files (./data/Focal_Tree_Neighborhood/). 

The subsequent section of code is focused on analyzing the how neighborhood trait similarities influence focal tree growth and survival. This section of code can be run using only the data that is provided in this git repository.
