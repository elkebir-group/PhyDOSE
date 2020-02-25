###########################################################################################
#
#     PhyDOSE version 1.0
#     Created by: Leah Weber, Nuraini Aguse, Nicholas Chia and  Mohammed El-Kebir
#
###########################################################################################
# 
# path_to_trees <- "/home/leah/Documents/UIUC/Research/SCS/PhyDOSE/PhyDOSE/data/patient2/trees_freqs"
# df_folder <-"/home/leah/Documents/UIUC/Research/SCS/PhyDOSE/PhyDOSE/data/patient2/distFeats"
# fmultiplier <- 2


###########################################################################################
#
#  Read in command line arguments
#
###########################################################################################
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 2){
  stop("At least one argument must be supplied (input file)", call.=FALSE)

}else{
  path_to_trees<- args[1]
  df_folder <- args[2]
}


home_path <- dirname(path_to_trees)


#read in optional arguments gamma and false negative rate
if(length(args) == 2){
  conf_level <- 0.95
  fn_rate <- 0
  fmult <- 1
}else if(length(args) == 3){
  conf_level <- 0.95
  fn_rate <- 0
  fmult <- args[3]
}else if(length(args) ==4){
  conf_level <- args[4]
  fn_rate <- 0
  fmult <- args[3]


}else{
  conf_level <- args[4]
  fn_rate <- args[5]
  fmult <- args[3]
}


conf_level <- as.numeric(conf_level)
fn_rate <- as.numeric(fn_rate)
fmult <- as.numeric(fmult)


  

###########################################################################################
#
#   Load dependies and helper scripts
#
###########################################################################################
library(pmultinom)
library(stringr)
library(dplyr,warn.conflicts=FALSE)
library(tidyr)


source('multinomial_union.r')
source('process_files.r')



###########################################################################################
#
#  Read in all trees in the input file and calculate the corresponding clonal prevalences
#
###########################################################################################

allTrees <- list.files(path_to_trees)

dist_feat_files <- list.files(df_folder)
phyDOSE <- data.frame(stringsAsFactors = F)
print("Calculating k*....")

for(a in allTrees){
  
  input <- main(file.path( path_to_trees,a), fmult)
  trees <- input$trees
  f_matrix <- input$f

  
  if(is.matrix(f_matrix)){
    
    samples <- nrow(f_matrix)
  }else{
    samples <- 1
  }
  
  u_list <- input$u
  
  tree_num <- as.numeric(str_extract(a, "[0-9]+")) 
  dist_feat_file <- dist_feat_files[str_detect(dist_feat_files, paste0("T", tree_num, "_"))]
  
  for(i in 1:samples){
    k_t <- calc_cells(df_path = file.path(df_folder, dist_feat_file), 
                      sample = i, 
                      tree_num = 1, 
                      u_list = u_list,
                      fn_rate = fn_rate,
                      gamma  = conf_level)
    
    tree_df <- data.frame(tree = a, tree_num = tree_num, sample = i, cells = k_t, stringsAsFactors = F)
    phyDOSE <- bind_rows(tree_df, phyDOSE)
  }
  
  
}







###########################################################################################
#
#   Determine k* by calculating the maximum per sample over all trees and then selecting
#    the sample with the minimum k*
#
###########################################################################################

phyDOSE.OUT <- phyDOSE %>% group_by(sample) %>%
  summarize(cells = max(cells)) %>%
  ungroup() %>%
  summarize(k_star = min(cells)) %>% 
  ungroup()

k_star <- phyDOSE.OUT$k_star


print(paste0("With confidence level ", conf_level, " and false negative rate ", fn_rate, ", k* = " , k_star))



       
     