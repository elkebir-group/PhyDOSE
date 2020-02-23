#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) <= 2){
  stop("At least two argument must be supplied (input file)", call.=FALSE)
  
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



#takes in dataframe wth all the trees  scs_file_name, number of cells to draw)
#returns a boolean vector with whether or not that tree was identified by that sample

scs_exp <- function(scs.dat,df_folder_path, name_map){
  

    

    cells <- ceiling(cells)
    #num_success <- rep(0, nrow(tree_df))
    scs.all <- data.frame()
    tree_support <- data.frame()
    for(i in 1:trials){
  
      scs.samples <- sample_n(scs.dat, cells, replace = bool_replace)
      
      for(j in 1:nrow(tree_df)){
        tree <- tree_df[j,]
        
        dfFamily <- df_family(file.path(df_folder_path, tree$tree_file_name))
        list.df <- convert_dfFamily_to_scs(dfFamily, colnames(scs.samples))
       pres <- unlist(lapply(list.df, is_present, scs.samples))
       
       ############Calculate primary support#######################3
       supported_df <- data.frame() 
       for(k in 1:length(pres)){
          if(pres[k] > 0){
            supported_df <- bind_rows(list.df[[k]], supported_df)
          }
       }
       supported_df <- supported_df %>% distinct()
        if(nrow(supported_df) > 0){
          support.prim <- count_cells(supported_df,scs.samples )
        }else{
          support.prim <- 0
        }
        
       ##########Calculate secondary support################
       support.sec <- numeric(length(dfFamily)) 
       for(k in 1:length(pres)){
         df <- list.df[[k]]
         
         feats_seen <- numeric(nrow(df))
         if(pres[k] > 0){
           for(l in 1:nrow(df)){
             feats_seen[l] <- nrow(semi_join(scs.samples, df[l,], by= colnames(scs.samples)))
           }
           support.sec[k] <- min(feats_seen)
           
         }
         
         
       }
      
     
        df <-unlist(lapply(dfFamily, paste, collapse = ":"))
        df_res <- as.data.frame(df) 
        df_res$df_id =1:length(df)
        df_res$success <- pres
        df_res$support.prime <- support.prim
        df_res$support.sec <- support.sec
        df_res$trial <- i
        df_res$tree_file_name <- unique(tree$tree_file_name)
        df_res$tree_num <- unique(tree$tree_num)
        df_res$cells <- cells
      
        #print(tail(df_res))
        tree_support <- bind_rows(df_res, tree_support)
        
        
      }
   
      
      
      scs.all <- bind_rows(scs.samples %>% mutate(trial = i), scs.all)
    }
    

    
    
    
    

  scs.all$stat <- stat
  scs.all$fn <- fn

 write.csv(scs.all,write_file, row.names=F)
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
