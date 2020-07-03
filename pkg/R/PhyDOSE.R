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
        return(list(k_star=-1, graphs=NULL))
    }
    phyDOSE <- data.frame(stringsAsFactors = F)
    resu_all <- main(filename, fmult)
    resu <- resu_all$main
    graphs <- resu_all$graphs
    #Rgraphviz::plot(graphs[[1]])
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
                              tree_num = i, 
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
    
    phyDOSE.FILTER <- (phyDOSE[phyDOSE$cells == k_star][1,])
    best_sample <- phyDOSE.FILTER$sample
    #print(phyDOSE.FILTER)
    #phyDOSE.TREE <- (phyDOSE[phyDOSE$cells == k_star][1,])$tree_num
    best_tree <- phyDOSE.FILTER$tree_num
    print(best_tree)
    
    # generate data for plot
    gamma_list <- seq(0.0, 0.99, by=0.001)
    cells <- numeric(length(gamma_list))
    index <- 1
    curr_df <- dff[[1]]
    input <- resu[[1]]
    u_list <- input$u
    for(g in gamma_list){
        cell_list <- calc_cells(curr_df, 
                                sample = best_sample, 
                                tree_num = best_tree, 
                                u_list = u_list,
                                fn_rate = fn_rate,
                                gamma  = g)
        cells[index] = cell_list
        
        index = index + 1
    }
    res <- data.frame(k=cells, gamma = gamma_list)
    
    
    # return all relevant information
    ret_all <- list(k_star=k_star, graphs=graphs, plot_vals=res, dff=dff)
    return (ret_all)
}