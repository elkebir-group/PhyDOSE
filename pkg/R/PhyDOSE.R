####################################################################
#
#   PhyDOSE parses filename and returns the number k* of cells to sequence
#
#   @param filename - Name of file containing trees and frequency matrix 
#                     using format specified in __
#   @param conf_level - Confidence level of successfully identifying the 
#                       true tree. Defaults to 0.95.
#   @param fn_rate - Expected false negative rate in single-cell sequencing
#                    experiments. Defaults to 0.
#   @param fmult - 
#
####################################################################

PhyDOSE <- function(filename, conf_level=0.95, fn_rate=0, fmult=1){
    library(pmultinom)
    library(stringr)
    library(dplyr,warn.conflicts=FALSE)
    library(tidyr)
    enum_err = FALSE
    dff <- tryCatch({
        enumerateDF(filename)
        }, 
             error = function(err){
                 print(paste("MY_ERROR:  ",err))
                 return(list())}
             )
    if(length(dff) == 0){
        return(-1)
    }
    phyDOSE <- data.frame(stringsAsFactors = F)
    resu <- main(filename, fmult)
    for(i in 1:length(resu)){
        input <- resu[[i]]
        curr_df <- dff[[i]]
        trees <- input$trees
        f_matrix <- input$f
        if(is.matrix(f_matrix)){
            
            samples <- nrow(f_matrix)
        }else{
            samples <- 1
        }
        u_list <- input$u
        for(j in 1:samples){
            k_t <- calc_cells(curr_df, 
                              sample = j, 
                              tree_num = 1, 
                              u_list = u_list,
                              fn_rate = fn_rate,
                              gamma  = conf_level)
            tree_name <- paste("tree",i)
            tree_df <- data.frame(tree = tree_name, 
                                  tree_num = i, 
                                  sample = j, 
                                  cells = k_t, 
                                  stringsAsFactors = F)
            phyDOSE <- bind_rows(tree_df, phyDOSE)
        }
    }
    phyDOSE.OUT <- phyDOSE %>% group_by(sample) %>%
        summarize(cells = max(cells)) %>%
        ungroup() %>%
        summarize(k_star = min(cells)) %>% 
        ungroup()
    
    k_star <- phyDOSE.OUT$k_star
    return (k_star)
}