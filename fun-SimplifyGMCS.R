# 
# # debuging
# gc()
# cat("\f")
# rm(list = ls())
# setwd("C:/Users/lvalcarcel/OneDrive - Tecnun/Postdoc/002 - MM Human-GEM NY - CoMMpass - CCLE EBI 2021-06/00 - GENERATE APP/gMCStool-shinyApp-v5/CODE_TEST")
# 
# # source("../helpers.R")
# 
# library(Matrix)
# 
# 
# gMCS.info <- new.env()
# # load("Data/gMCSs_EssentialTasks_CultureBiomass_combined_HumanGEMv1.4.1_ENSEMBL.rdata", envir = gMCS.info)
# load("../Data/gMCSs_all_cases_HumanGEMv1.4.0_ENSEMBL.rdata", envir = gMCS.info)
# # gMCS.info <- gMCS.info$gMCS.info.raw$EssentialTasks_CultureMedium
# 
# # gMCS.info$gMCSs.ENSEMBL.txt.list <- lapply(strsplit(gMCS.info$gMCSs.ENSEMBL.txt, '--'),sort) # array of unique gMCSs along task, but separated
# # gMCS.info$gMCSs.ENSEMBL.txt.list.num <- lapply(lapply(strsplit(gMCS.info$gMCSs.ENSEMBL.txt, '--'),sort),match,gMCS.info$genes.gMCSs.ENSEMBL) # array of unique gMCSs along task, but separated
# 
# 
# # idx <- rowSums(gMCS.info$gMCSs.ENSEMBL.mat) <= 6
# # gMCS.info$gMCSs.ENSEMBL.mat <- gMCS.info$gMCSs.ENSEMBL.mat[idx,]        # matrix that relates gMCS and genes
# # gMCS.info$genes.gMCSs.ENSEMBL <- gMCS.info$genes.gMCSs.ENSEMBL    # all genes contained in the gMCS
# # gMCS.info$gMCSs.ENSEMBL.txt <- gMCS.info$gMCSs.ENSEMBL.txt[idx]        # array of unique gMCSs along tasks, text. genes separated by "--"
# # gMCS.info$table.genes.HumanGEM <- gMCS.info$table.genes.HumanGEM  # table that relates genes in ENSEMBL, SYMBOL and ENTREZ_ID
# # gMCS.info$gMCSs.ENSEMBL.txt.list <- gMCS.info$gMCSs.ENSEMBL.txt.list[idx]
# # gMCS.info$gMCSs.ENSEMBL.txt.list.num <- gMCS.info$gMCSs.ENSEMBL.txt.list.num[idx]
# 
# 
# 
# 
# 
# # trace = T
# # gMCS.info <- gMCS.info$gMCS.info.raw$EssentialTasks_CultureMedium


FunSimplifyGMCS <- function(gMCS.info,  # gMCS information
                            complete = T, # not only remove essential genes, perform complete analysis
                            trace = T){  #if true, information is printed during the running of step.
  
  # browser()
  # extract info from the gMCSs
  gMCSs.ENSEMBL.mat <- gMCS.info$gMCSs.ENSEMBL.mat        # matrix that relates gMCS and genes
  genes.gMCSs.ENSEMBL <- gMCS.info$genes.gMCSs.ENSEMBL    # all genes contained in the gMCS
  gMCSs.ENSEMBL.txt <- gMCS.info$gMCSs.ENSEMBL.txt        # array of unique gMCSs along tasks, text. genes separated by "--"
  if ("gMCSs.ENSEMBL.txt.list" %in% names(gMCS.info)){
    gMCSs.ENSEMBL.txt.list <- gMCS.info$gMCSs.ENSEMBL.txt.list
    gMCSs.ENSEMBL.txt.list.num <- gMCS.info$gMCSs.ENSEMBL.txt.list.num
  } 
  table.genes.HumanGEM <- gMCS.info$table.genes.HumanGEM  # table that relates genes in ENSEMBL, SYMBOL and ENTREZ_ID
  
  
  
  ###########################  simplify the list   #############################
  
  # we can smplify the list by eliminating those gMCSs that contain gMCSs of 
  # lower order. Which implies that they are not minimal
  
  if (trace){
    print(paste0("Dimmensions of gMCSs.ENSEMBL.mat are (",paste(dim(gMCSs.ENSEMBL.mat),collapse = ", "),")"))
  }
  
  gMCSs.ENSEMBL.mat.raw <- gMCSs.ENSEMBL.mat
  rownames(gMCSs.ENSEMBL.mat) <- 1:nrow(gMCSs.ENSEMBL.mat)
  
  # order by length
  gMCSs.ENSEMBL.length <- rowSums(gMCSs.ENSEMBL.mat)
  gMCSs.ENSEMBL.mat <- gMCSs.ENSEMBL.mat[order(gMCSs.ENSEMBL.length, decreasing = FALSE),]
  
  
  # eliminate gMCS that contains gMCS of order 1
  gMCSs.ENSEMBL.length <- rowSums(gMCSs.ENSEMBL.mat)
  gMCSs.ENSEMBL.mat <- gMCSs.ENSEMBL.mat[order(gMCSs.ENSEMBL.length, decreasing = FALSE),]
  
  idx_length_1 <- which(gMCSs.ENSEMBL.length==1)
  gmcs_order_1 <- colSums(gMCSs.ENSEMBL.mat[idx_length_1, , drop = F])
  gmcs_order_1 <- names(gmcs_order_1)[gmcs_order_1>0]
  idx_length_not_1 <- which(gMCSs.ENSEMBL.length > 1)
  
  aux_mat = gMCSs.ENSEMBL.mat[idx_length_not_1, gmcs_order_1, drop = F] # auxiliary matrix of the genes of these gMCS
  idx = idx_length_not_1[rowSums(aux_mat)>0] # if sum is 1 or more, this is because the gmcs contains an essential gene
  
  if (length(idx)>0){ 
    gMCSs.ENSEMBL.mat <- gMCSs.ENSEMBL.mat[-idx,,drop = F]
  }
  
  
  
  # reduce by removing thoe gmcs of order (k) that contains gmcs of order (1,...,k-1)
  # to do so, we search for order k, gMCSs of order k+1,...,inf that contains 
  # the gmcs of order k
  
  if (complete) {
    genes.gMCSs.ENSEMBL <- colnames(gMCSs.ENSEMBL.mat) # reduced list of genes in gMCSs
    gMCSs.ENSEMBL.length <- rowSums(gMCSs.ENSEMBL.mat) # calculate new lenght of gMCSs
    
    k = 1
    while( k < nrow(gMCSs.ENSEMBL.mat)){
      # print(k)
      gmcs <- which(gMCSs.ENSEMBL.mat[k,]>0) # columns that contain genes of that gMCS
      
      idx <- which(gMCSs.ENSEMBL.length>length(gmcs))
      if (length(idx)>0){ 
        aux_mat <- gMCSs.ENSEMBL.mat[idx, gmcs, drop = F] # auxiliary matrix of the genes of these gMCS
        idx2 <- rowMeans(aux_mat) # how many genes are in each gmcs that are part of gmcs_k
        idx <- idx[idx2==1] # if mean is 1, this is because the gmcs contains gmcs_k
        
        if (length(idx)>0){ 
          gMCSs.ENSEMBL.mat <- gMCSs.ENSEMBL.mat[-idx,] # remove the rows
          gMCSs.ENSEMBL.mat <- gMCSs.ENSEMBL.mat[,colSums(gMCSs.ENSEMBL.mat)>=1] # remove empty columns
          genes.gMCSs.ENSEMBL <- colnames(gMCSs.ENSEMBL.mat) # reduced list of genes in gMCSs
          gMCSs.ENSEMBL.length <- rowSums(gMCSs.ENSEMBL.mat) # calculate new lenght of gMCSs
        }
      }
      k = k+1 # increase k
      # print(paste0("Dimmensions of gMCSs.ENSEMBL.mat have been reduced to (",paste(dim(gMCSs.ENSEMBL.mat),collapse = ", "),")"))
    }
  }
  
  # Transform to text
  gMCSs.ENSEMBL.txt <- apply(gMCSs.ENSEMBL.mat,1, function(x){paste0(genes.gMCSs.ENSEMBL[x>0.5], collapse = "--")})
  gMCSs.ENSEMBL.txt <- as.character(gMCSs.ENSEMBL.txt)
  
  # browser()
  genes.gMCSs.SYMBOL <- table.genes.HumanGEM$SYMBOL[match(genes.gMCSs.ENSEMBL, table.genes.HumanGEM$ENSEMBL)]
  gMCSs.SYMBOL.txt <- apply(gMCSs.ENSEMBL.mat,1, function(x){paste0(genes.gMCSs.SYMBOL[x>0.5], collapse = "--")})
  gMCSs.SYMBOL.txt <- as.character(gMCSs.SYMBOL.txt)
  
  genes.gMCSs.ENSEMBL <- colnames(gMCSs.ENSEMBL.mat) # reduced list of genes in gMCSs
  gMCSs.ENSEMBL.length <- rowSums(gMCSs.ENSEMBL.mat) # calculate new lenght of gMCSs
  
  if (trace){
    print(paste0("Dimmensions of gMCSs.ENSEMBL.mat have been reduced to (",paste(dim(gMCSs.ENSEMBL.mat),collapse = ", "),")"))
  }
  
  
  
  
  ###########################   Export the results   ###########################
  # browser()
  
  if (exists("gMCSs.ENSEMBL.txt.list") & setequal(colnames(gMCSs.ENSEMBL.mat.raw),colnames(gMCSs.ENSEMBL.mat))) {
    gMCSs.ENSEMBL.txt.list <- gMCSs.ENSEMBL.txt.list[as.numeric(as.character(rownames(gMCSs.ENSEMBL.mat)))]
    gMCSs.ENSEMBL.txt.list.num <- gMCSs.ENSEMBL.txt.list.num[as.numeric(as.character(rownames(gMCSs.ENSEMBL.mat)))]
  } else {
    gMCSs.ENSEMBL.txt.list <- lapply(strsplit(gMCSs.ENSEMBL.txt, '--'),sort) # array of unique gMCSs along task, but separated,
    gMCSs.ENSEMBL.txt.list.num <- lapply(gMCSs.ENSEMBL.txt.list,match,genes.gMCSs.ENSEMBL) # array of unique gMCSs along task, but separated,
    
  }
  
  return(list(
    gMCSs.ENSEMBL.mat = gMCSs.ENSEMBL.mat,
    genes.gMCSs.ENSEMBL = genes.gMCSs.ENSEMBL,
    gMCSs.ENSEMBL.length = gMCSs.ENSEMBL.length,
    gMCSs.ENSEMBL.txt = gMCSs.ENSEMBL.txt,
    gMCSs.SYMBOL.txt = gMCSs.SYMBOL.txt,
    # gMCSs.ENSEMBL.txt.list = gMCSs.ENSEMBL.txt.list,
    gMCSs.ENSEMBL.txt.list.num = gMCSs.ENSEMBL.txt.list.num, 
    
    # previous information
    table.genes.HumanGEM = gMCS.info$table.genes.HumanGEM,
    table.gMCSs = gMCS.info$table.gMCSs,
    gMCSs.ENSEMBL = gMCS.info$gMCSs.ENSEMBL,
    fullname = gMCS.info$fullname
  ))
  
  
  
  
}


# for (n in names(gMCS.info$gMCS.info.raw)[1:4]){
#   START <- Sys.time()
#   print(paste(n))
#   
#   FunSimplifyGMCS(gMCS.info$gMCS.info.raw[[n]])
#   
#   print(paste(n, Sys.time()-START))
# }

