#!/usr/bin/env Rscript


#Rscript PhyDOSE.r [path to SCS csv file][path to distinguishing features directory][path to mutation name mapping file]
args = commandArgs(trailingOnly=TRUE)



if(length(args) < 2){
  stop("At least two argument must be supplied (single-cell data and path to distinguishing features)", call.=FALSE)
  
}else{
  scsFile <- args[1]
  df_folder <- args[2]
  
}

if(length(args) > 2){
  name_changes <-TRUE
  nameMapFile <- args[3]
}else{
  name_changes <- FALSE
}



library(dplyr,warn.conflicts=FALSE)
library(stringr)
source('process_files.r')

scs.dat <- read.csv(scsFile, stringsAsFactors = F)


#if different names are used for the mutations in the SCS experiments than in the original
# tree files, then there must me a name map provided in csv format that maps the 
# scs column names (first column) to the mutation names used in the trees (second column)
if(name_changes){
  nameMap <- read.csv(nameMap, stringsAsFactors = F)
  col.df <- data.frame(scs_names =colnames(scs.dat), stringsAsFactors = F)
  col.df <- col.df %>% inner_join(nameMap)
  to_keep <- filter(col.df, tree_names!= "")
  scs.dat<- select(scs.dat, one_of(to_keep$scs_names))
  colnames(scs.dat) <- to_keep$tree_names
}





scs_exp <- function(scs.dat,df_folder_path){
   dist_feat_files <- list.files(df_folder_path)
  tree_support <- data.frame()
  for(f in dist_feat_files){
    dfFamily <- df_family(file.path(df_folder_path, f))
    list.df <- convert_dfFamily_to_scs(dfFamily, colnames(scs.dat))
    pres <- unlist(lapply(list.df, is_present, scs.dat))
    
    
    ############Calculate support #######################3
    supported_df <- data.frame() 
    for(k in 1:length(pres)){
      if(pres[k] > 0){
        supported_df <- bind_rows(list.df[[k]], supported_df)
      }
    }
    supported_df <- supported_df %>% distinct()
    if(nrow(supported_df) > 0){
      support <- count_cells(supported_df,scs.dat)
    }else{
      support <- 0
    }
    
    df_tree <- data.frame(tree = f, support = support, stringsAsFactors = F)
    tree_support  = bind_rows(df_tree, tree_support)
  }
  
  
  tree_support$tree_number <- as.numeric(str_replace(str_extract(tree_support$tree, "_T[0-9]*"), "_T", ""))
  tree_support <- select(tree_support, tree, tree_number, support) %>% arrange(-support) %>% filter(support > 0)
  
  return(tree_support)

}


convert_dfFamily_to_scs <- function(dfList, cnames){
  data.list <- lapply(dfList, convert_to_vector, cnames)
  return(data.list)
}

convert_to_vector <- function(vec, cnames){
  if(substr(cnames[1],0,1) == "X"){
    muts <- as.numeric(str_extract(cnames, '[0-9]+'))
  }
  else{
    muts <- cnames
  }
  #
  names.df <- data.frame(n =cnames)
  vec2 <- str_split(vec, pattern = " ", simplify = FALSE)
  for(v in vec2){
    x <- ifelse(muts %in% v, 1, 0)
    names.df <- cbind(names.df,x)
  }
  names.df <- as.data.frame(t(names.df[,-1])) 
  colnames(names.df) <- cnames
  rownames(names.df) <- NULL
  
  return(names.df)
  
  
}


count_cells <- function(df, scs.dat){
  df_support <- nrow(semi_join(scs.dat, df, by= colnames(scs.dat)))
  
  return(df_support)
  
}

is_present <- function(df, scs.dat){
  temp <- semi_join(df, scs.dat, by = colnames(scs.dat))
  pres <- ifelse(nrow(temp) == nrow(df), 1, 0)
  
  return(pres)
  
}



tree_support <- scs_exp(scs.dat, df_folder)

cat("Tree","\t", "Support", "\n")

for(i in 1:nrow(tree_support)){
  cat(tree_support[i,"tree_number"],"\t", tree_support[i,"support"], "\n" )
}

