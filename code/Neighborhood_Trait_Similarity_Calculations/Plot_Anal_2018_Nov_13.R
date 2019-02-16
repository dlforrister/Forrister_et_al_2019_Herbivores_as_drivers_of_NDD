# This scripts serves as the backbone all data analysis. It calls the following for scripts in order to 1) Load in required trait data, 2) Calculate growth between census intervals 3) load in plot analysis function designed to calculate neighborhood trait and herbivores similarities for all focal trees.

source("./code/Neighborhood_Trait_Similarity_Calculations/Analysis_Setup.R") 
source("./code/Neighborhood_Trait_Similarity_Calculations/plot_growth_dataframe.R")
source("./code/Neighborhood_Trait_Similarity_Calculations/plot_anal_function_2018_Nov_13.R")

cen_int_function <- function(start_cen,end_cen,cen_int, fixed_null = F) {
  plot_growth_df(start_cen,end_cen,cen_int)
  
 #Reg Plot analysis
  c4.ppp <- ppp(c4$gx,c4$gy,window=owin_plot)
  marks(c4.ppp) <- c4$treeID 
  
  #actual plot data
  plot_anal(plot_df = c4,ppp_object = c4.ppp, database_name = 'BCI_INGA_2018_Nov_13')
  
  #generation of random data for null model
  repeat{
  rand.ppp<-rpoispp(lambda = (sum_plot$intensity-0.03), max(quad_count), owin_plot)
  if(rand.ppp$n < length(c4$treeID)){ break }}
  marks(rand.ppp) <- sample(c4$treeID,rand.ppp$n) 
  rand_plot_df <- c4[which(c4$treeID %in% rand.ppp$marks),] # this makes it so only trees that were randomly     picked are listed in C4
  
  
  #null function. 
  #When the flag fixed_null is set to True (fixed_null == T) both a fixed null model, where only the species names are swapped among individuals in the plot, and a non fixed model where both the location and the species names are randomly generated
  
  plot_anal(plot_df = rand_plot_df,ppp_object = rand.ppp, database_name = "BCI_INGA_Null_2018_Nov_13")
  if(fixed_null == T) {
    ingas <- which(c4$sp %in% names$sp)
    rand_plot_df <- c4
    rand_plot_df$sp[which(rand_plot_df$sp %in% names$sp)] <- rand_plot_df$sp[ingas[shuffle(ingas, control = how(within = Within(type = "series")))]]
  plot_anal(plot_df = rand_plot_df,ppp_object = rand.ppp, database_name = "BCI_INGA_Null_fixed_2018_Nov_13")}
  }  

#The above function gets run for each cencus interval  
cen_int_function(bci.full1,bci.full2,1,T)
cen_int_function(bci.full2,bci.full3,2,T)
cen_int_function(bci.full3,bci.full4,3,T)
cen_int_function(bci.full4,bci.full5,4,T)
cen_int_function(bci.full5,bci.full6,5,T)
cen_int_function(bci.full6,bci.full7,6,T)
cen_int_function(bci.full7,bci.full8,7,T)