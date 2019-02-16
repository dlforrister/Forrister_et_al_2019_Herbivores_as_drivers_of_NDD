source("./code/Growth_Survival_Analysis/BCI_NULL_Growth_Survival_10_meter.R")
library(lmerTest)
Null_Growth_Survival(fixed_null = T)

BCI_sap <- subset(BCI_Norm,dbh <50)

par(mfrow = c(5,2))

plot(BCI_sap$cg_nondef_w ~ BCI_sap$gap_TL, main = round(summary(lm(BCI_sap$cg_nondef_w ~ BCI_sap$gap_TL))$coefficients[2,4],3)) |abline(lm(BCI_sap$cg_nondef_w ~ BCI_sap$gap_TL),col = "red") | 

plot(BCI_sap$cg_nondef_l ~ BCI_sap$gap_TL, main = round(summary(lm(BCI_sap$cg_nondef_l ~ BCI_sap$gap_TL))$coefficients[2,4],3)) | abline(lm(BCI_sap$cg_nondef_l ~ BCI_sap$gap_TL),col = "red")

plot(BCI_sap$cg_nondef_e ~ BCI_sap$gap_TL, main = round(summary(lm(BCI_sap$cg_nondef_e ~ BCI_sap$gap_TL))$coefficients[2,4],3)) | abline(lm(BCI_sap$cg_nondef_e ~ BCI_sap$gap_TL),col = "red")

plot(BCI_sap$cg_nondef_h ~ BCI_sap$gap_TL, main = round(summary(lm(BCI_sap$cg_nondef_h ~ BCI_sap$gap_TL))$coefficients[2,4],3)) |abline(lm(BCI_sap$cg_nondef_h ~ BCI_sap$gap_TL),col = "red")

