source("./code/Neighborhood_Trait_Similarity_Calculations/Inga_Trait_Measures_2.R")
source("./code/Neighborhood_Trait_Similarity_Calculations/herbivore_similarity_calculation.R")
library(picante)
library(tidyr)

### get quadrats:
load("./data/Tree_Tables/bci.tree1.rdata")
load("./data/Tree_Tables/bci.tree2.rdata")
load("./data/Tree_Tables/bci.tree3.rdata")
load("./data/Tree_Tables/bci.tree4.rdata")
load("./data/Tree_Tables/bci.tree5.rdata")
load("./data/Tree_Tables/bci.tree6.rdata")
load("./data/Tree_Tables/bci.tree7.rdata")
load("./data/Tree_Tables/bci.tree8.rdata")
census <- list(bci.tree1,bci.tree2,bci.tree3,bci.tree4,bci.tree5,bci.tree6,bci.tree7,bci.tree8)
c = 1
s = 1

web_stat_life_stage <- data.frame(NRI = NA,NRT.obs = NA,NRI.rand = NA, NTI = NA, NTI.obs=NA, NTI.rand=NA,life_stage= NA, census = NA)
web_stat_life_stage <- web_stat_life_stage[-1,]
for(c in 1:length(census)){

census_x = census[[c]]  
census_x$x_quad <- cut(as.numeric(census_x$gx),breaks = 20, labels = 1:20)
census_x$y_quad <- cut(as.numeric(census_x$gy),breaks = 10, labels = 1:10)
census_x$quadrat <- paste(census_x$x_quad,census_x$y_quad,sep="-")
census_x_a <- subset(census_x,status == "A")
census_x_a$gx <- jitter(census_x_a$gx)
census_x_a$gy <- jitter(census_x_a$gy)
census_x_a$size_class <- cut(census_x_a$dbh,c(9,50,100,4000), labels = c("sap","juv","adult"))
inga_sp <- c("ingaco","ingafa","ingago","ingama","ingape","ingaqu","ingas1","ingasa","ingaum") #turn this of ones included in analysis. 
census_x_a_inga <- census_x_a[which(census_x_a$sp %in% inga_sp),]
census_x_a_inga$sp <- as.factor(census_x_a_inga$sp)

quad_sp_count <- function(quad){
  quadrat <- subset(census_x_a_inga,quadrat == quad)
  samp_size <- table(quadrat$sp)
  return(samp_size)
  }

#picante NRI.
for(s in 1:length(unique(census_x_a_inga$size_class))){
    census_x_a_inga_sc <- subset(census_x_a_inga, size_class == unique(census_x_a_inga$size_class)[s])
    quadtrats <- as.matrix(t(sapply(unique(census_x_a_inga_sc$quadrat), FUN = quad_sp_count)))
    webb_stat_NRI <- ses.mpd(quadtrats,herb_dist_jaccard_leps_saw_ish_norm, null.model = "taxa.labels",abundance.weighted = T,iterations = 999)
    webb_stat_NTI <- ses.mntd(quadtrats,herb_dist_jaccard_leps_saw_ish_norm, null.model = "taxa.labels",abundance.weighted = T,iterations = 999)
    temp <- data.frame(NRI = webb_stat_NRI$mpd.obs.z, NRI.obs = webb_stat_NRI$mpd.obs, NRI.rand =webb_stat_NRI$mpd.rand.mean, NTI = webb_stat_NTI$mntd.obs.z, NTI.obs = webb_stat_NTI$mntd.obs, NTI.rand =webb_stat_NTI$mntd.rand.mean,life_stage = unique(census_x_a_inga$size_class)[s], census = c(1:8)[c])
    web_stat_life_stage <- rbind(web_stat_life_stage,temp)
  }
}

sap <- subset(web_stat_life_stage, life_stage == "sap")
sap.t.test.NRI <- t.test(x = sap$NRI.obs, y= sap$NRI.rand)
sap.t.test.NTI <- t.test(x = sap$NTI.obs, y= sap$NTI.rand)

juv <- subset(web_stat_life_stage, life_stage == "juv")
juv.t.test.NRI <- t.test(x = juv$NRI.obs, y= juv$NRI.rand)
juv.t.test.NTI <- t.test(x = juv$NTI.obs, y= juv$NTI.rand)

adult <- subset(web_stat_life_stage, life_stage == "adult")
adult.t.test.NRI <- t.test(x = adult$NRI.obs, y= adult$NRI.rand)
adult.t.test.NTI <- t.test(x = adult$NTI.obs, y= adult$NTI.rand)

length(unique(bci.tree1$treeID))


long_DF <- web_stat_life_stage %>% gather(Measure, Value, NRI:NTI)

head(long_DF)
boxplot(Value ~ Measure*life_stage,data = long_DF)
abline(h = 0, col = "red")

 



