source("./code/Growth_Survival_Analysis/BCI_NULL_Growth_Survival_10_meter.R")
library(lmerTest)
Null_Growth_Survival(fixed_null = T)

BCI_sap <- subset(BCI_Norm,dbh <50)
head(BCI_sap)
dep_vars <- names(BCI_sap)[c(13:25)]
names_dep_vars <- c( "BA:conspecific density","BB:congeneric density", "Def:Chemistry","Def:Ants","Def:Development","Def:Hairs","Def:Phenology_pearson","CA:Resource:Wood Density","CB:Resource:Elem.Comp","CC:Resource:Height","CD: Resource:Leaf_Morph","E:Herb_Sim_bray_leps_saw_ish","E:Herb_Sim_jaccard_leps_saw_ish","A: gap_specialist_TAK_PDC")

length(dep_vars)
length(names_dep_vars)
shading <- data.frame(min = seq(from = 0.5, to = max(as.numeric(as.factor(names_dep_vars))), by = 1),
                      max = seq(from = 1.5, to = max(as.numeric(as.factor(names_dep_vars))) + 0.5, by = 1),
                      col = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1))

#### Combined model accounting for gap:
results_growth_gap = data.frame(NA,NA,NA,NA,NA,NA)
results_growth_gap <- results_growth_gap[-1,]
names(results_growth_gap) <- c("Estimate"," Std. Error","pval","max","min","sig")

head(growth)
for( n in dep_vars) {
  growth <- BCI_sap[-which(is.na(BCI_sap$norm_grth)),]
  library(splitstackshape)
  growth <- as.data.frame(cSplit(growth, "id_cen", sep = "_"))
  dep_var <- growth[,n]
  
  #the particular model was changed between the below options to generate the data generated in the supplementary  material
  
  #model <- lm(growth$norm_grth ~ dep_var) #Does not take into account number of gap specialists (Table S5)
  #model<-lmer(growth$norm_grth ~ dep_var + (1|growth$gap_TL) + (1|growth$quadrat) + (1|growth$id_cen_1)) #Model that accounts for spatial autocorrelation by including quadrat ID and census interval as random effects (Table S8)
  model<-lmer(growth$norm_grth ~ dep_var + (1|growth$gap_TL)) #Model that was used to generate Figure 1.
  row <- as.data.frame(summary(model)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
  confinf_int <-  confint(model, "dep_var", level = 0.95)
  row$max = confinf_int[2]
  row$min = confinf_int[1]
  row$sig = cut(row$`Pr(>|t|)`,breaks = c(-0.001,0.05,1),labels = c("sig","non_sig"))
  results_growth_gap <-rbind(results_growth_gap,row["dep_var",])
}

model<-lm(norm_grth ~ gap_TL,data=BCI_sap)
row <- as.data.frame(summary(model)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
confinf_int <-  confint(model, "gap_TL", level = 0.95)
row$max = confinf_int[2]
row$min = confinf_int[1]
row$sig = cut(row$`Pr(>|t|)`,breaks = c(-0.001,0.05,1),labels = c("sig","non_sig"))
results_growth_gap <-rbind(results_growth_gap,row["gap_TL",])

results_growth_gap$var <- names_dep_vars

#adjust for multiple inference
results_growth_gap$adjust_p <- p.adjust(results_growth_gap$`Pr(>|t|)`,method ="holm")
results_growth_gap$adjust_p_sig <- cut(results_growth_gap$adjust_p,breaks = c(-0.001,0.05,1),labels = c("sig","non_sig"))


ggplot(results_growth_gap, aes(x = var, y = Estimate, ymin = min, ymax = max,shape=adjust_p_sig)) + theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=20, face='bold',hjust = 0.5),
        axis.title=element_text(size=18,face="bold"),
        axis.text=element_text(size=18),
        legend.key.size =  unit(0.3, "in"),
        axis.text.x  = element_text(angle=90, vjust=0.5),
        legend.text = element_text(size=18))+ 
  geom_errorbar(colour="black", position=position_dodge(width=0.75)) +
  geom_point(position=position_dodge(width=0.75),size = 5) + 
  geom_hline(yintercept = 0, col = "red") +
  labs(title= "Effect of Neighborhood Similarity on Focal Tree Growth", x=  "Neighborhood Similarity Type", y= "Effect Size (sd)") + 
  geom_rect(inherit.aes=FALSE,data = shading, aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, fill = factor(col)), alpha = 0.25) + 
  scale_fill_manual(values = c("transparent", "gray53"), guide = FALSE) + scale_shape_manual(breaks=c("sig","marg_sig","non_sig"),values=c(15,16,17))


#### Survival model with Gap

results_surv_gap = data.frame(NA,NA,NA,NA,NA,NA,NA)
results_surv_gap <- results_surv_gap[-1,]
names(results_surv_gap) <- c("Estimate"," Std. Error","pval","max","min","sig","model_type")

BCI_sap$Survival <- 0
BCI_sap$Survival[which(BCI_sap$life_event == "S")] <- 1
BCI_sap <- as.data.frame(cSplit(BCI_sap, "id_cen", sep = "_"))
summary(model)
head(BCI_sap)
for( n in dep_vars) {
  dep_var <- BCI_sap[,n]
  #model<- glm(BCI_sap$Survival ~ dep_var,family = "binomial")
  #model<- glmer(BCI_sap$Survival ~ dep_var + (1|BCI_sap$quadrat) + (1|BCI_sap$id_cen_1),family = "binomial")
  model<- glmer(BCI_sap$Survival ~ dep_var + (1|BCI_sap$gap_TL),family = "binomial" )
  row <- as.data.frame(summary(model)$coefficients[,c("Estimate","Std. Error","Pr(>|z|)")])
  confinf_int <-  confint(model, "dep_var", level = 0.95)
  row$max = confinf_int[2]
  row$min = confinf_int[1]
  row$sig = cut(row$`Pr(>|z|)`,breaks = c(-0.001,0.05,1),labels = c("sig","non_sig"))
  results_surv_gap <-rbind(results_surv_gap,row["dep_var",])
}

model<-glm(Survival ~  gap_TL ,data=BCI_sap,family = "binomial" )
row <- as.data.frame(summary(model)$coefficients[,c("Estimate","Std. Error","Pr(>|z|)")])
confinf_int <-  confint(model, "gap_TL", level = 0.95)
row$max = confinf_int[2]
row$min = confinf_int[1]
row$sig = cut(row$`Pr(>|z|)`,breaks = c(-0.001,0.05,1),labels = c("sig","non_sig"))
results_surv_gap <-rbind(results_surv_gap,row["gap_TL",])

results_surv_gap$var <- names_dep_vars

results_surv_gap$adjust_p <- p.adjust(results_surv_gap$`Pr(>|z|)`,method ="holm")
results_surv_gap$adjust_p_sig <- cut(results_surv_gap$adjust_p,breaks = c(-0.001,0.05,1),labels = c("sig","non_sig"))


ggplot(results_surv_gap, aes(x = var, y = Estimate, ymin = min, ymax = max,shape=adjust_p_sig)) + theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=20, face='bold',hjust = 0.5),
        axis.title=element_text(size=18,face="bold"),
        axis.text=element_text(size=18),
        legend.key.size =  unit(0.3, "in"),
        axis.text.x  = element_text(angle=90, vjust=0.5),
        legend.text = element_text(size=18))+ 
  geom_errorbar(colour="black", position=position_dodge(width=0.75)) +
  geom_point(position=position_dodge(width=0.75),size = 5) + 
  geom_hline(yintercept = 0, col = "red") +
  labs(title= "Effect of Neighborhood Similarity on Focal Tree Survival \n Model controling for Gap specialist", x=  "Neighborhood Similarity Type", y= "Effect Size (sd)") + 
  geom_rect(inherit.aes=FALSE,data = shading, aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, fill = factor(col)), alpha = 0.25) + 
  scale_fill_manual(values = c("transparent", "gray53"), guide = FALSE) + scale_y_continuous(limits = c(-0.25, 0.25))


