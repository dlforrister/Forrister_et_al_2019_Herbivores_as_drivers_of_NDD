plot_anal <- function(plot_df,ppp_object,database_name) {
  inga_focal <- which(plot_df$sp %in% names$sp) #list of trees of the genus inga... list of rows in c2 that are inga
  inga_focal <- (as.data.frame(ppp_object[which(ppp_object$marks %in% plot_df$treeID[inga_focal]),])) 
  foc_loc <- inga_focal[,c(3,1,2)]
  names(foc_loc) <- c("treeID","gx","gy")
  
  foc_loc_L <-foc_loc[which(foc_loc$gx >=12.6),] # all focal trees greater than 12.6 meters from L edge
  foc_loc_R <-foc_loc_L[which(foc_loc_L$gx <= 987.4),] # all focal trees less than 12.6 meters from R edge
  foc_loc_T <-foc_loc_R[which(foc_loc_R$gy <=487.4 ),] # all focal trees greater than 12.6 meters from top edge
  foc_loc_B <-(foc_loc_T[which(foc_loc_T$gy >= 12.6 ),]) # all focal trees greater than 12.6 meters from top edge
  
  library(doParallel)
  library(foreach)
  
  cores=detectCores()
  cl <- makeCluster(cores[1])
  clusterEvalQ(cl, {
    library(RMySQL)
    mydb = dbConnect(MySQL(), user='u6000251', password='********', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
    NULL
  })
  registerDoParallel(cl)
  #
all_ftsummary <-foreach(i=1:nrow(foc_loc_B),.noexport = c("mydb"), .export=c("An","lab","fields","chem_sim_sedio_norm","ants_sim_norm","dev_sim_norm","hair_sim_norm", "phen_pearson_sim_norm","wood_sim_norm","element_sim_norm" ,"height_sim_norm","leaf_sim_norm","TL_gap_sp", "herb_sim_norm_jaccard_leps_saw_ish"),.errorhandling=c("pass") ,.packages = c("spatstat","RMySQL")) %dopar% {
    ft_id <- foc_loc_B[i,1]
    pot <- spatstat::disc(radius= 12.6, centre=c(foc_loc_B[i,2],foc_loc_B[i,3]))
    focal <-ppp_object[pot,] 
    pair_dist <- pairdist(focal, squared = T)
    row.names(pair_dist) <- focal$marks
    col <- which(focal$marks == ft_id)
    dist_matrix <- as.data.frame(pair_dist[,col])
    names(dist_matrix) <- "r2"
    dist_matrix$r <- sqrt(dist_matrix$r2)
    dist_matrix$pair_ID <- focal$marks
    dist_matrix$dbh <- plot_df$dbh[which(plot_df$treeID %in% focal$marks)]
    dist_matrix$basal_area <- plot_df$basal_area[which(plot_df$treeID %in% focal$marks)]
    dist_matrix$sp <- plot_df$sp[which(plot_df$treeID %in% focal$marks)]
    dist_matrix$an <- cut(dist_matrix$r, breaks = An,labels = lab[order(lab)],right = F)
    dist_matrix$size_class <-cut(dist_matrix$dbh, c(9,30,50,70,90,110,130,150,170,190,3000), labels= 1:10)
    dist_matrix$sap <- cut(dist_matrix$dbh, c(9,50,4000),labels= c(1,0))
    dist_matrix$sap <- as.numeric(levels(dist_matrix$sap)[dist_matrix$sap])
    dist_matrix$con_spec <- ifelse(grepl(plot_df$sp[which(plot_df$treeID == ft_id)],dist_matrix$sp),1,0)
    dist_matrix$con_gen <- ifelse(grepl("inga",dist_matrix$sp),1,0)
    dist_matrix$con_gen <- ifelse(grepl("ingami",dist_matrix$sp),0,dist_matrix$con_gen)
    dist_matrix$con_gen[grepl(plot_df$sp[which(plot_df$treeID == ft_id)],dist_matrix$sp)] <- 0
    dist_matrix$het <- ifelse(dist_matrix$con_spec + dist_matrix$con_gen != 0,0,1)
    dist_matrix$mort <- F
    dist_matrix$mort[which(plot_df$life_event[which(plot_df$treeID %in% focal$marks)] == "D")] <- T
    dist_matrix$rec <- F
    dist_matrix$rec[which(plot_df$life_event[which(plot_df$treeID %in% focal$marks)] == "R")] <- T
    dist_matrix$norm_growth <- plot_df$norm_grth[which(plot_df$treeID %in% focal$marks)]
    dist_matrix$norm_agb <- plot_df$norm_agb[which(plot_df$treeID %in% focal$marks)]
    dist_matrix <-dist_matrix[-which(dist_matrix$pair_ID == ft_id),] # when we've computed all the variables, one shoudl remove the neighbor of distance 0 because that is the distance to itself
   
    #defense
    for(n in 1:length(dist_matrix$sp))
      if(length(which(names(chem_sim_sedio_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_chem_s[n] <- chem_sim_sedio_norm[,which(names(chem_sim_sedio_norm) == dist_matrix$sp[n])][which(row.names(chem_sim_sedio_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_chem_s[n] <- 0}
    
    for(n in 1:length(dist_matrix$sp))
      if(length(which(names(ants_sim_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_ants[n] <- ants_sim_norm[,which(names(ants_sim_norm) == dist_matrix$sp[n])][which(row.names(ants_sim_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_ants[n] <- 0}
    
    for(n in 1:length(dist_matrix$sp))
      if(length(which(names(dev_sim_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_dev[n] <- dev_sim_norm[,which(names(dev_sim_norm) == dist_matrix$sp[n])][which(row.names(dev_sim_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_dev[n] <- 0}
    
    for(n in 1:length(dist_matrix$sp))
      if(length(which(names(hair_sim_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_hair[n] <- hair_sim_norm[,which(names(hair_sim_norm) == dist_matrix$sp[n])][which(row.names(hair_sim_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_hair[n] <- 0}

  #phenology #2 circular stats peason
    
    for(n in 1:length(dist_matrix$sp))
    if(length(which(names(phen_pearson_sim_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_phen_pearson[n] <- phen_pearson_sim_norm[,which(names(phen_pearson_sim_norm) == dist_matrix$sp[n])][which(row.names(phen_pearson_sim_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_phen_pearson[n] <- 0}
    
#resource
    
    for(n in 1:length(dist_matrix$sp))
      if(length(which(names(wood_sim_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_wood[n] <- wood_sim_norm[,which(names(wood_sim_norm) == dist_matrix$sp[n])][which(row.names(wood_sim_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_wood[n] <- 0}
    
    for(n in 1:length(dist_matrix$sp))
      if(length(which(names(element_sim_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_element[n] <- element_sim_norm[,which(names(element_sim_norm) == dist_matrix$sp[n])][which(row.names(element_sim_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_element[n] <- 0}
    
    for(n in 1:length(dist_matrix$sp))
      if(length(which(names(height_sim_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_height[n] <- height_sim_norm[,which(names(height_sim_norm) == dist_matrix$sp[n])][which(row.names(height_sim_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_height[n] <- 0}
    
    for( n in 1:length(dist_matrix$sp))
      if(length(which(names(leaf_sim_norm) == dist_matrix$sp[n])) > 0) {dist_matrix$congen_leaf[n] <- leaf_sim_norm[,which(names(leaf_sim_norm) == dist_matrix$sp[n])][which(row.names(leaf_sim_norm) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$congen_leaf[n] <- 0}
    
    #gap
    for( n in 1:length(dist_matrix$sp))
      if(length(which(TL_gap_sp == dist_matrix$sp[n])) > 0 & dist_matrix$dbh[n] < 50) {dist_matrix$gap_TL[n] <- 1} else{dist_matrix$gap_TL[n] <- 0}
    
 #herb_comparisons
    #leps_saw_include_single_host_jaccard
    for( n in 1:length(dist_matrix$sp))
    if(length(which(names(herb_sim_norm_jaccard_leps_saw_ish) == dist_matrix$sp[n])) > 0) {dist_matrix$herb_sim_norm_jaccard_leps_saw_ish[n] <- herb_sim_norm_jaccard_leps_saw_ish[,which(names(herb_sim_norm_jaccard_leps_saw_ish) == dist_matrix$sp[n])][which(row.names(herb_sim_norm_jaccard_leps_saw_ish) == plot_df$sp[which(plot_df$treeID == ft_id)])]} else{dist_matrix$herb_sim_norm_jaccard_leps_saw_ish[n] <- 0}
    
    ftsummary <-as.data.frame(ft_id)
    ftsummary$iteration <- i
    ftsummary$census <- plot_df$census[which(plot_df$treeID == ft_id)]
    ftsummary$sp <- plot_df$sp[which(plot_df$treeID == ft_id)]
    ftsummary$quadrat <- plot_df$quadrat[which(plot_df$treeID == ft_id)]
    ftsummary$gx <- plot_df$gx[which(plot_df$treeID == ft_id)]
    ftsummary$gy <- plot_df$gy[which(plot_df$treeID == ft_id)]
    ftsummary$life_event <- plot_df$life_event[which(plot_df$treeID == ft_id)]
    ftsummary$dbh <- plot_df$dbh[which(plot_df$treeID == ft_id)]
    ftsummary$basal_area <- plot_df$basal_area[which(plot_df$treeID == ft_id)]
    ftsummary$sc <- plot_df$size_class[which(plot_df$treeID == ft_id)]
    ftsummary$norm_agb <- plot_df$norm_agb[which(plot_df$treeID == ft_id)]
    ftsummary$norm_grth <- plot_df$norm_grth[which(plot_df$treeID == ft_id)]
    
    ftsummary <-do.call("rbind", replicate(20, ftsummary, simplify = FALSE))
    ftsummary$an <- lab[order(lab)]
    ftsummary$tot <- aggregate(dist_matrix$basal_area, by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cs <-aggregate((dist_matrix$basal_area*dist_matrix$con_spec), by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cg <- aggregate((dist_matrix$basal_area*dist_matrix$con_gen), by=list(Category=dist_matrix$an), FUN=sum)$x
   
    ftsummary$cg_def_c_s <-aggregate((dist_matrix$basal_area*dist_matrix$congen_chem_s), by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cg_def_a <-aggregate((dist_matrix$basal_area*dist_matrix$congen_ants), by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cg_def_d <-aggregate((dist_matrix$basal_area*dist_matrix$congen_dev), by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cg_def_h <-aggregate((dist_matrix$basal_area*dist_matrix$congen_hair), by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cg_def_p_pearson <-aggregate((dist_matrix$basal_area*dist_matrix$congen_phen_pearson), by=list(Category=dist_matrix$an), FUN=sum)$x

    ftsummary$cg_nondef_w <-aggregate((dist_matrix$basal_area*dist_matrix$congen_wood), by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cg_nondef_e <-aggregate((dist_matrix$basal_area*dist_matrix$congen_element), by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cg_nondef_h <-aggregate((dist_matrix$basal_area*dist_matrix$congen_height), by=list(Category=dist_matrix$an), FUN=sum)$x
    ftsummary$cg_nondef_l <-aggregate((dist_matrix$basal_area*dist_matrix$congen_leaf), by=list(Category=dist_matrix$an), FUN=sum)$x
    
    ftsummary$cg_herb_jaccard_leps_saw_ish <-aggregate((dist_matrix$basal_area*dist_matrix$herb_sim_norm_jaccard_leps_saw_ish), by=list(Category=dist_matrix$an), FUN=sum)$x
    
    ftsummary$gap_TL <- (aggregate((dist_matrix$gap_TL), by=list(Category=dist_matrix$an), FUN=sum)$x/aggregate((dist_matrix$gap_TL), by=list(Category=dist_matrix$an), FUN=length)$x)*100
    
    dbWriteTable(mydb,database_name,ftsummary, field.types=fields,row.names=F,overwrite=F, append=T)
    print(i)
  }
  stopCluster(cl)
}

