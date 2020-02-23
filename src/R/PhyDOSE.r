###########################################################################################
#
#     PhyDOSE version 1.0
#     Created by: Leah Weber, Nuraini Aguse, Nicholas Chia and  Mohammed El-Kebir
#
###########################################################################################


###########################################################################################
#
#  Read in command line arguments
#
###########################################################################################
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  stop("At least one argument must be supplied (input file)", call.=FALSE)

}else{
  treeFile <- args[1]
}


home_path <- dirname(treeFile)
tree_file_name <- basename(treeFile)

#read in optional arguments gamma and false negative rate
if(length(args) == 1){
  conf_level <- 0.95
  fn_rate <- 0
  df_folder <- file.path(home_path, "distFeats")
}else if(length(args) == 2){
  conf_level <- args[2]
  fn_rate <- 0
  df_folder <- file.path(home_path, "distFeats")
}else if(length(args) ==3){
  conf_level <- args[2]
  fn_rate <- args[3]
  df_folder <- file.path(home_path, "distFeats")

}else{
  conf_level <- args[2]
  fn_rate <- args[3]
  df_folder <- args[4]
}


conf_level <- as.numeric(conf_level)
fn_rate <- as.numeric(fn_rate)


  

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
input <- main(treeFile)

trees <- input$trees
f_matrix <- input$f

if(is.matrix(f_matrix)){

  samples <- nrow(f_matrix)
}else{
  samples <- 1
}

u_list <- input$u

dist_feat_files <- list.files(df_folder)

phyDOSE <- expand.grid(tree =dist_feat_files, sample= 1:samples)
phyDOSE$tree_num <-  as.numeric(str_replace(str_extract(phyDOSE$tree, "_T[0-9]*"), "_T", ""))


###########################################################################################
#
#   Run PhyDOSE for each tree and sample 
#
###########################################################################################



phyDOSE <- phyDOSE %>% rowwise() %>% mutate(cells = calc_cells(df_path = file.path(df_folder, tree), 
                                                    sample = sample, 
                                                    tree_num = tree_num, 
                                                    u_list = u_list,
                                                    fn_rate = fn_rate,
                                                    gamma  = conf_level)) %>% 
  mutate(gamma = conf_level, fn = fn_rate, tree_set = tree_file_name) %>%
  select(tree_set, gamma, fn, tree_num, sample, cells) %>% ungroup()


###########################################################################################
#
#   Determine k* by calculating the maximum per sample over all trees and then selecting
#    the sample with the minimum k*
#
###########################################################################################

phyDOSE.OUT <- phyDOSE %>% group_by(tree_set, gamma, fn, sample) %>%
  summarize(cells = max(cells)) %>%
  group_by(tree_set, gamma, fn) %>%
  summarize(k_star = min(cells)) %>% 
  ungroup()

k_star <- phyDOSE.OUT$k_star


print(paste0("With confidence level ", conf_level, " and false negative rate ", fn_rate, ", k* = " , k_star))


       
     