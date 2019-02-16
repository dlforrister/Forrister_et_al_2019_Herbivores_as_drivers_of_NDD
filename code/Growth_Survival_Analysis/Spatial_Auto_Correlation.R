source("./code/Growth_Survival_Analysis/BCI_NULL_Growth_Survival_10_meter.R")
library(lmerTest)
Null_Growth_Survival(fixed_null = T)

BCI_sap <- subset(BCI_Norm,dbh <50)
growth <- BCI_sap[-which(is.na(BCI_sap$norm_grth)),]

library(rgdal)
library(spdep)


#convert dataframe into spatial dataframe
coordinates(growth) <- growth[,c("gx","gy")]

coords <- coordinates(growth)

#calculate spatial weights matrix
growth_1_nb = knn2nb(knearneigh(coords, k = 5))
dists = unlist(nbdists(growth_1_nb, coords))
growth_1_nb_2 = dnearneigh(coords, d1 = 0, d2 = 20)
growth_1_lw_B = nb2listw(growth_1_nb_2, style='B')

#test for spatial autocorrelation in model.

#model_1: original used in manuscript

model<-lmer(growth$norm_grth ~ dep_var + (1|growth$gap_TL))
moran.mc(residuals(model), growth_1_lw_B, nsim=999)

#model_2: added quadrat and census interval as random effect.

model_quad<-lmer(growth$norm_grth ~ dep_var + (1|growth$gap_TL) + (1|growth$quadrat) + (1|growth$))
moran.mc(residuals(model_quad), growth_1_lw_B, nsim=999)


#spatial lag models for everything

#These take a very long time to run.


growth_lm_splag_cs = lagsarlm(growth$norm_grth ~ growth$cs + (1|growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cs)

growth_lm_splag_cg = lagsarlm(growth$norm_grth ~ growth$cg + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg)

growth_lm_splag_cg_def_c_s = lagsarlm(growth$norm_grth ~ growth$cg_def_c_s + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_def_c_s)

growth_lm_splag_cg_def_a = lagsarlm(growth$norm_grth ~ growth$cg_def_a + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_def_a)

growth_lm_splag_cg_def_d = lagsarlm(growth$norm_grth ~ growth$cg_def_d + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_def_d)

growth_lm_splag_cg_def_h = lagsarlm(growth$norm_grth ~ growth$cg_def_h + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_def_h)

growth_lm_splag_cg_def_p_pearson = lagsarlm(growth$norm_grth ~ growth$cg_def_p_pearson + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_def_p_pearson)

growth_lm_splag_cg_nondef_w = lagsarlm(growth$norm_grth ~ growth$cg_nondef_w + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_nondef_w)

growth_lm_splag_cg_nondef_e = lagsarlm(growth$norm_grth ~ growth$cg_nondef_e + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_nondef_e)

growth_lm_splag_cg_nondef_h = lagsarlm(growth$norm_grth ~ growth$cg_nondef_h + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_nondef_h)

growth_lm_splag_cg_nondef_l = lagsarlm(growth$norm_grth ~ growth$cg_nondef_l + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_nondef_l)

growth_lm_splag_cg_herb_jaccard_leps_saw_ish = lagsarlm(growth$norm_grth ~ growth$cg_herb_jaccard_leps_saw_ish + growth$gap_TL, growth, growth_1_lw_B)
summary(growth_lm_splag_cg_herb_jaccard_leps_saw_ish)



spat_model = data.frame(model = c("cs","chem","ants","dev","hairs","pearson","wooddensity","elements","hight","leaf","herbivores"),pval = c(0.008362,0.006013,0.2932,0.2886,0.8262,0.2892,0.4777,0.4647,0.312,0.3518,0.01038))
spat_model$p_adjust <- p.adjust(spat_model$pval,method ="holm")

