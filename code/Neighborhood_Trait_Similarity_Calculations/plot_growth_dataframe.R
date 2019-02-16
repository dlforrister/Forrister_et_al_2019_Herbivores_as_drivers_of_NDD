#!/usr/bin/Rscript
plot_growth_df <- function(start_cen,end_cen,cen_int) {
  c1<-start_cen[complete.cases(start_cen[,c("treeID","gx","gy")]),] # all trees with tree ID and Xy cordinates in 
  cens_2 <- end_cen[complete.cases(end_cen[,c("treeID","gx","gy")]),c("treeID","status","dbh")]
  c2 <- merge(cens_2,c1,by = "treeID", all.y =T)
  c2$life_event <-paste(c2$status.y,c2$status.x,sep="_") #status.y = census 1
  c2$dbh.x[which(is.na(c2$dbh.x))] <- 0
  c2$dbh.y[which(is.na(c2$dbh.y))] <- 0
  c2$dbh <- apply(c2[,c("dbh.x","dbh.y")], 1, max)
  
  if(length(which(c2$life_event == "A_AD")) >0) { c2 <- c2[-which(c2$life_event == "A_AD"),]}
  if(length(which(c2$life_event == "A_AM")) >0){ c2 <- c2[-which(c2$life_event == "A_AM"),] }
  if(length(which(c2$life_event == "AD_A")) >0){ c2 <- c2[-which(c2$life_event == "AD_A"),]}
  if(length(which(c2$life_event == "AD_AD")) >0){ c2 <- c2[-which(c2$life_event == "AD_AD"),]}                            
  if(length(which(c2$life_event == "AM_A")) >0){ c2 <- c2[-which(c2$life_event == "AM_A"),]} 
  if(length(which(c2$life_event == "AR_AR")) >0){ c2 <- c2[-which(c2$life_event == "AR_AR"),]}
  if(length(which(c2$life_event == "D_D")) >0){ c2 <- c2[-which(c2$life_event == "D_D"),]} 
  if(length(which(c2$life_event == "P_P")) >0){ c2 <- c2[-which(c2$life_event == "P_P"),]}
  if(length(which(c2$life_event == "AR_A")) >0){ c2 <- c2[-which(c2$life_event == "AR_A"),]}
  if(length(which(c2$life_event == "A_M")) >0){ c2 <- c2[-which(c2$life_event == "A_M"),]}
  if(length(c2$life_event[which(c2$life_event == "M_D")]) > 1) {c2$life_event[which(c2$life_event == "M_D")] <- "A_D"}
  
  c2$life_event <- factor(c2$life_event,levels = c("A_A","A_D","P_A"),labels = c("S","D","R"))                
  c2$gx <- jitter(x = c2$gx)
  c2$gy <- jitter(x=c2$gy)
  
  owin_plot <<- owin(xrange=c(-1,1001), yrange=c(-1,501))
  c2.ppp <- ppp(c2$gx,c2$gy,window=owin_plot) #these are all points in the plot. Only P(precensus) has bene removed. Note, row names are kept, so row.name = treeID. 
  
  grate=individual_grow.table(cnsdata=list(start_cen,end_cen),mingrow=0.0,mindbh=10,maxdbh=10000,center=1992)
  grate$size_class <-cut(grate$dbh1, c(9,30,50,70,90,110,130,150,170,190,3000), labels= 1:10)
  grate$sp_sc <- paste(grate$species,"_",grate$size_class,sep="")
  
  grate$mean_growth <- ave(grate$CRGrowth,grate$sp_sc)
  grate$sd_growth <- ave(grate$CRGrowth,grate$sp_sc, FUN = sd)
  grate$norm_grth <- ((grate$CRGrowth-grate$mean_growth)/grate$sd_growth)
  
  c3 <- merge(c2,grate, by= "treeID",all.x = T)
  c3$mean_agb <- ave(c3$agb,c3$sp,FUN = mean)
  c3$sd_agb <- ave(c3$agb,c3$sp,FUN = sd)
  c3$norm_agb <- (c3$agb-c3$mean_agb)/c3$sd_agb
  c3$dbh[which(is.na(c3$dbh))] <- 0
  c3 <- c3[-which(c3$dbh == 0),]
  c3$size_class <-cut(c3$dbh, c(-1,9,30,50,70,90,110,130,150,170,190,3000), labels= 0:10)
  c3$census <- cen_int
  c3$basal_area <- round((pi*(c3$dbh/2)^2),0)
  c4 <- c3[,c("treeID","sp","quadrat","gx.x","gy.x","dbh","size_class","life_event","norm_agb","norm_grth","census","basal_area")]
  names(c4) <- c("treeID","sp","quadrat","gx","gy","dbh","size_class","life_event","norm_agb","norm_grth","census","basal_area")
  c4 <<- c4
}
