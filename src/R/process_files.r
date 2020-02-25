


##########################################################
#
#     Function to read in an F matrix from a standard file
#     and return a dataframe in long format matching how
#     the datafile is formmatted
#
#      @param charVec - a character vector of lines read in from a file
#
#      @param start - which line in the character vector should the file 
#                     start being read. Defaults to 1. 
# 
#
#
#
###########################################################
get_F_matrix <- function(charVec, start=1){
  fmatrix <- data.frame()
  fmat <- F
  for(i in 1:length(charVec)){
    if(str_detect(charVec[i], "sample_index")){
      
      charVec[i] <- str_replace(charVec[i],"#", "")
      
      header <- str_split(charVec[i], "\t")
      header <- header[[1]]
      header <- c(header[1:6], c("fminus", "fplus"))
      
      fmat <- T
      
      
      next
    }
    if(fmat){
      new_vec <- str_replace_all(charVec[i], ";","-")
      new_vec <- str_replace_all(charVec[i], "\n","")
      new_vec <- str_split(new_vec, "\t")
      new_vec <- data.frame(val =new_vec[[1]], stringsAsFactors = F)
      if(length(new_vec) !=0){
        new_vec$header <- header
        new_vec <- new_vec %>% pivot_wider( names_from = "header", values_from="val")
        #print(new_vec)
        fmatrix <- rbind(fmatrix, new_vec)
      }
      
      
      
      
    }
  }
  
  return(fmatrix)
}


#################################################################################
# 
#   A helper function to determine all of the ancesters for all nodes
#
#
#   @param vecs - a vector of the names of all internal vertices
#
#   @param root_node - specifies which character label is the root node
#
#    returns a list of all the ancesters for each node
#
#
##################################################################################
create_parent_list <- function(vecs, root_node){
  
  parents_all <- list()
 
  branches <- str_split(vecs, pattern = " ")
  for(b in branches){
    
    for(i in 1:length(b)){
      
      if(!(b[i] %in% names(parents_all))){
        parents_all[b[i]] <- root_node
      }
      
      #it's a child so add the immediate parent
      if( i > 1 && b[1] != root_node){
        
        parents_all[[b[i]]] <- c(parents_all[[b[i]]], b[1])
      }
    }
  }
 

continue <- T  
while(continue){
    cur_len <- lapply(parents_all, length)
    for(p in names(parents_all)){

      if( p != root_node){
        for(b in parents_all[[p]]){
          
          if(b != root_node){
            parents_all[[p]] <- unique(c(parents_all[[p]], parents_all[[b]]))
          }
        }
      
        
      }
      
      
    }
    
    new_len <- lapply(parents_all, length)

    if(identical(new_len, cur_len)){
      break;
    }
    
  }
  
parents_all[[root_node]] <- NA
  
  return(parents_all)
  
}






###############################################################################
#
#
#    Creates a tree in binary matrix form by reading in specified lines of a file 
#      containing one or more trees in standard DOT form
#
#   @param lines - character vector of all lines read from a file
#
#   @param start - the line at which the tree starts in the file
#
#  @param end    - the line at which the tree ends in the file
#   
#   @param key  -  a dataframe mapping the character index to the character label
#
#     
#
#
################################################################################

create_tree <- function(lines, start, end, key){
  vecs <- rev(lines[start:end])

  treeMat <-  matrix(rep(0, nrow(key)**2), 
                     nrow=nrow(key), ncol=nrow(key), dimnames = list(clones= as.character(0:(nrow(key)-1)),mutations=key$character_index) )

  find_root <- key %>% arrange(desc(fplus)) %>% top_n(1, fplus)
  root_node <- find_root$character_label[1]
  
  root_node_number <- key$character_index[match(root_node, key$character_label)]

  treeMat[,root_node_number] <- 1
  #print(root_node)
  
  parent_list <- create_parent_list(vecs, root_node) 
  #print(parent_list)
  
  
  clones <- list()
  edge <- 1
  clones[[edge]] <- c(root_node)
  
  for(v in names(parent_list)){
    if(v == root_node){
      next
    }
    edge <- edge + 1
    new_clone <- v
    b <- key$character_index[which(key$character_label==v)]
    treeMat[edge,b] <-1
    for(p in parent_list[[v]]){
      new_clone <- c(new_clone,p)
      b <- key$character_index[which(key$character_label==p)]
      treeMat[edge,b] <-1
    }
    
    clones[[edge]] <- new_clone
    
  }
  clone_names <- vector()
  i <- 1
  for(c in clones){
    c <- sort(c)
    if(length(c) > 1){
      clone_names[i] <- paste(c, collapse = " ")
    }
    else{
      clone_names[i] <- c
    }
    i <- i+ 1
    
  }
  
  subkey <- select(key, -fplus) %>% distinct()
  #print(clone_names)
  #print(subkey$character_label)
  dimnames(treeMat) = list(clones= clone_names,mutations=subkey$character_label)
  return(treeMat)
}

###############################################################################
#
#
#    A wrapper around create tree that scans a file and finds where all
#     of the trees are located and calls create_tree on those lines
#
#   @param lines - character vector of all lines read from a file
#
#   
#   @param key  -  a dataframe mapping the character index to the character label
#
#
################################################################################s
get_trees <- function(lines, key){
  trees <- list()
  start_tree <- FALSE
  
  i <- 1
  tree_index <- 0
  while(!str_detect(lines[i], pattern= "characters")){
    
    
    if(str_detect(lines[i], pattern = "edges") &&  !start_tree){
      
      s <- i + 1
      start_tree <- TRUE
      i<- i  + 1
      next
      
    }
    
    if((str_detect(lines[i], pattern = "edges") || str_detect(lines[i], pattern = "sites")) 
       && start_tree){
      e <- i -1
      tree_index <- tree_index + 1
     # print(paste0("build tree ", s, " ", e))

    
      trees[[tree_index]] <- create_tree(lines, start = s ,end = e, key)
      s <- i + 1
      
    }
    
    
    
    i <- i + 1
    
    
  }
  
  
  return(trees)
}



############################################################
#   Given a tree matrix, B, and a frequency matrix F, calculate
#    the corresponding U matrix by multiply F times B inverse
#
#  @param treeMat - a binary matrix representing a tree
#
#  @param f - a VAF frequency matrix
#
#
#
#######################################################
calc_u <- function(treeMat, f){
  t_inv <- solve(treeMat)
  u <- f %*% t_inv
  return(u)
}



#####################################################
#
#   Helper function to convert the Frequency data in dataframe format
#    to a matrix
#   
#  @param df - a dataframe representing an F matrix
#
#  @param key - a dataframe mapping the character index to the character label
#
###################################################
clean_fmatrix <-function(df, key){
  fmatrix <- df %>% select(sample_index, character_label, fplus) %>%
    mutate(fplus = as.numeric(fplus)) %>%
    pivot_wider(id_cols = sample_index, names_from = character_label, values_from = fplus)  %>% 
    select(-sample_index)
  fmatrix <- as.matrix(fmatrix)
  col.order <- key$character_label
  fmatrix <- fmatrix[,col.order]
  fmatrix <- fmatrix
  return(fmatrix)
}


####################################################
#
#   Reads in a distinguishing feature family file and 
#    and returns the distinguishing feature family in list format
#
#    @param fname - the name including the path of the file
#
###################################################
df_family <- function(fname){
  lines <- readLines(fname)
  lines <- str_remove_all(lines, "'")
  df_list <- str_split(lines, pattern = ",")
  
  return(df_list)
  
}


########################################################
#
#  A helper function to rename the distinguishing features
#   to what they are named in the U matrix for easy lookup
#
#   @param  df_fam - a distinguishing feature family in list format
#
#   @param u - a U matrix where the colnames are the clone names
#
#
########################################################
rename_df <- function(df_fam, u){
  clones <- str_split(colnames(u), pattern = " ")
  df_rename <- list()
  
  for(i in 1:length(df_fam)){
    new_df_vector <- vector()
    df <- df_fam[[i]]
    for(featurette in df){
      
      featurette <- str_split(featurette, pattern = " ")[[1]]
      
      for(c in clones){
        
        m <- match(featurette, c)
        m <- m[!is.na(m)]
        if(length(m) == length(featurette)){
          #print("true")
          new_df_vector <-unique(c(new_df_vector,paste(c, collapse = " ")))
          
        }
      }
    }
    
    df_rename[[i]] <- new_df_vector
    
  }
  return(df_rename)
}






##############################################
#  Main funtion to process a tree file and return a list containing
#   all pertinent matrices (f,u and b)
#
#  @param fname - the name of the tree and f file to be processed
#  @param path - the path where the file is located. Defaults to the empty string
#  @paramm adjustClusters - currently not in use but could be used to adjust for false negatives
#
#
##################################################
main <- function(fname, fmultiplier=1 ){

  fil <- readLines(fname)
  fmat <- get_F_matrix(fil)
  key <- filter(fmat, sample_index ==0) %>% 
    select( character_index, character_label, fplus) %>% distinct() 
  fmatrix <- clean_fmatrix(fmat, key)
  
  if(!is.matrix(fmatrix)){
    fmatrix <- t(as.matrix(fmatrix))
  }


  trees <- get_trees(fil, key)
         
  key <- key %>% arrange(desc(fplus))
  col.order <- key$character_label
  fmatrix <- fmatrix[, col.order]
  for(i in 1:length(trees)){
    t <- trees[[i]] 
    t <- t[, col.order]
    trees[[i]] <- t
  }
  fmatrix <- fmatrix * fmultiplier
  u_list <- lapply(trees, calc_u, f= fmatrix)
  res<- list(name = fname,  f= fmatrix, u = u_list, trees = trees)
  
  return(res)
}



