CalculateEssentialGenes_singleTH <- function(gene.exp, gMCS.info, sample.class, singleTH, nWorkers, threshold_logFC = 1e-2){
  
  #browser()
  
  # prepare cluster
  library(doSNOW)
  cl <- makeCluster(nWorkers)
  registerDoSNOW(cl)
  
  
  # extract info
  gMCSs.ENSEMBL.mat <- gMCS.info$gMCSs.ENSEMBL.mat        # matrix that relates gMCS and genes
  genes.gMCSs.ENSEMBL <- gMCS.info$genes.gMCSs.ENSEMBL    # all genes contained in the gMCS
  gMCSs.ENSEMBL <- gMCS.info$gMCSs.ENSEMBL                # list of tasks that contains lists of gMCS
  gMCSs.ENSEMBL.txt <- gMCS.info$gMCSs.ENSEMBL.txt        # array of unique gMCSs along tasks, text. genes separated by "--"
  table.genes.HumanGEM <- gMCS.info$table.genes.HumanGEM  # table that relates genes in ENSEMBL, SYMBOL and ENTREZ_ID
  
  # Calculate first and second most expressed gene in each gMCS in each sample ####
  gene.exp[is.na(gene.exp)] <- 0 # eliminate nans
  
  # dim(gene.exp)
  
  nn <- dim(gMCSs.ENSEMBL.mat)[1]
  
  gene.first.exp = gene.first.ENSEMBL =
    gene.second.exp = gene.second.ENSEMBL =
    gene.first.exp.2 = gene.first.ENSEMBL =
    gene.second.exp.2 = gene.second.ENSEMBL =
    matrix(NA,nrow = nn, ncol = dim(gene.exp)[2])
  
  colnames(gene.first.exp) = colnames(gene.first.ENSEMBL) = 
    colnames(gene.second.exp) = colnames(gene.second.ENSEMBL) = 
    colnames(gene.exp)
  
  index_gmcs_length_1 <- as.numeric(which(rowSums(gMCSs.ENSEMBL.mat)==1))
  
  k <- list()
  
  gene.exp.gMCSs <- gene.exp[genes.gMCSs.ENSEMBL,]
  rownames(gene.exp.gMCSs) <- genes.gMCSs.ENSEMBL
  gene.exp.gMCSs[is.na(gene.exp.gMCSs)] <- 0
  aux_zeros_mat <- matrix(0, nrow = 3, ncol = dim(gene.exp)[2])
  
  print(paste("Processing the",nn, "gMCSs for",ncol(gene.exp), "samples:"))
  
  
  obtain_1_2_3_genes_gMCS <- function(i = 1){
    # i = 15000
    # i = 2
    # gmcs <- gMCSs.ENSEMBL.list[[i]]
    gmcs <- unlist(strsplit(gMCSs.ENSEMBL.txt[i], '--'))
    if (length(gmcs)==1){
      gMCSs.ENSEMBL.exp <- rbind(matrix(gene.exp[gmcs,], nrow = 1), aux_zeros_mat)
      gMCSs.ENSEMBL.exp <- apply(gMCSs.ENSEMBL.exp,2,unlist)
      colnames(gMCSs.ENSEMBL.exp) <- colnames(gene.exp)
    } else {
      gMCSs.ENSEMBL.exp <- rbind(as.matrix(gene.exp[gmcs,]), aux_zeros_mat)
    }
    gmcs <- c(gmcs, rep("", nrow(aux_zeros_mat)))
    rownames(gMCSs.ENSEMBL.exp) <- gmcs
    
    # save expresion of most expressed genes
    aux <- apply(gMCSs.ENSEMBL.exp,2,function(x){sort(x, decreasing = T)[1:2]})
    # save which gene is the most expressed
    aux2 <- apply(gMCSs.ENSEMBL.exp,2,function(x){order(x, decreasing = T)[1:2]})
    
    return(list(first.exp = aux[1,],
                second.exp = aux[2,],
                first.ENSEMBL = gmcs[aux2[1,]],
                second.ENSEMBL = gmcs[aux2[2,]]))
  }
  
  
  pb <- txtProgressBar(max = nn, style = 3)
  progress <- function(i) {setTxtProgressBar(pb, i); incProgress(1/nn)}
  opts <- list(progress = progress)
  withProgress(message = 'Calculating most expressed gene by gMCS and sample',
               detail = '\nThis may take a while...', {
                 results <- foreach(i=1:nn, .options.snow = opts) %dopar% obtain_1_2_3_genes_gMCS(i) 
               })
  
  gene.first.exp <- do.call(rbind,lapply(results, function(x){matrix(x$first.exp, nrow = 1)}))
  gene.second.exp <- do.call(rbind,lapply(results, function(x){matrix(x$second.exp, nrow = 1)}))
  gene.first.ENSEMBL <- do.call(rbind,lapply(results, function(x){matrix(x$first.ENSEMBL, nrow = 1)}))
  gene.second.ENSEMBL <- do.call(rbind,lapply(results, function(x){matrix(x$second.ENSEMBL, nrow = 1)}))
  
  colnames(gene.first.exp) = colnames(gene.first.ENSEMBL) = 
    colnames(gene.second.exp) = colnames(gene.second.ENSEMBL) = 
    colnames(gene.exp)
  
  
  for (sample in 1:dim(gene.exp)[2]) {
    idx_not_expressed <- which(gene.first.exp[,sample]<1e-3)
    gene.first.exp[idx_not_expressed,sample] <- 0
    gene.first.ENSEMBL[idx_not_expressed,sample] <- ""
    
    idx_not_expressed <- which(gene.second.exp[,sample]<1e-3)
    gene.second.exp[idx_not_expressed,sample] <- 0
    gene.second.ENSEMBL[idx_not_expressed,sample] <- ""
    
    
    # filter 
    gene.second.exp[index_gmcs_length_1, sample] <- 0
    gene.second.ENSEMBL[index_gmcs_length_1, sample] <- ""
  }
  # })
  
  gene.second.exp[index_gmcs_length_1,] <- 0
  gene.second.ENSEMBL[index_gmcs_length_1,] <- ""
  
  
  
  # generate ratio over single threshold
  gene.first.ratio <- gene.first.exp / singleTH
  gene.second.ratio <- gene.second.exp / singleTH
  gene.ratio <- gene.exp / singleTH
  
  
  # Calculate essential genes by sample ####
  gene.first.log2ratio <- log2(gene.first.ratio)
  gene.second.log2ratio <- log2(gene.second.ratio)
  
  
  genes.ENSEMBL.list <- list()
  for( n in names(gMCSs.ENSEMBL)){
    genes.ENSEMBL.list[[n]] <- unique(as.character(as.matrix(gMCSs.ENSEMBL[[n]])))
    genes.ENSEMBL.list[[n]] <- genes.ENSEMBL.list[[n]][genes.ENSEMBL.list[[n]]!=""]
    genes.ENSEMBL.list[[n]] <- genes.ENSEMBL.list[[n]][!is.na(genes.ENSEMBL.list[[n]])]
    genes.ENSEMBL.list[[n]]
  }
  
  
  genes.ENSEMBL.essential.list <- list()
  for( n in names(gMCSs.ENSEMBL)) {
    genes.ENSEMBL.essential.list[[n]] <- gMCSs.ENSEMBL[[n]][apply(gMCSs.ENSEMBL[[n]],1,function(x){sum(x!="")})==1,]
    genes.ENSEMBL.essential.list[[n]] <- unique(as.character(as.matrix(genes.ENSEMBL.essential.list[[n]])))
    genes.ENSEMBL.essential.list[[n]] <- genes.ENSEMBL.essential.list[[n]][genes.ENSEMBL.essential.list[[n]]!=""]
    genes.ENSEMBL.essential.list[[n]] <- genes.ENSEMBL.essential.list[[n]][!is.na(genes.ENSEMBL.essential.list[[n]])]
    genes.ENSEMBL.essential.list[[n]]
  }
  
  unlist(lapply(genes.ENSEMBL.list, length))
  unlist(lapply(genes.ENSEMBL.essential.list, length))
  
  
  genes.ENSEMBL.essential <- gMCSs.ENSEMBL.txt[!grepl('--', gMCSs.ENSEMBL.txt)]
  
  
  # if the gMCS is only one gene, make it essential     ####
  gene.first.log2ratio[rowSums(gMCSs.ENSEMBL.mat)==1,] <- (+1)
  gene.second.log2ratio[rowSums(gMCSs.ENSEMBL.mat)==1,] <- (-1)
  
  
  ################################################################
  ####    single essential genes                              ####
  ################################################################
  
  
  genes.gMCSs.ENSEMBL.12 <- unique(c(unique(gene.first.ENSEMBL), unique(gene.second.ENSEMBL)))
  genes.gMCSs.ENSEMBL.12 <- genes.gMCSs.ENSEMBL.12[genes.gMCSs.ENSEMBL.12!=""]
  length(genes.gMCSs.ENSEMBL.12)
  
  essentail.single.gMCS <- gene.first.log2ratio > (+threshold_logFC) & gene.second.log2ratio < (-threshold_logFC)
  
  num.essential.gene <- as.data.frame(matrix(0,nrow = length(genes.gMCSs.ENSEMBL),ncol = length(levels(sample.class))+2))
  colnames(num.essential.gene) <- c("gen","num.gMCS",levels(sample.class))
  rownames(num.essential.gene) <- genes.gMCSs.ENSEMBL
  num.essential.gene$gen <- genes.gMCSs.ENSEMBL
  
  mat.essential.gene <- matrix(0,nrow = length(genes.gMCSs.ENSEMBL),ncol = dim(gene.exp)[2])
  colnames(mat.essential.gene) <- colnames(gene.exp)
  rownames(mat.essential.gene) <- genes.gMCSs.ENSEMBL
  
  
  print(sum(essentail.single.gMCS))
  
  # define function to obtain essential genes ####
  
  obtain_essential_single_genes_gMCS <- function(i){
    gen = genes.gMCSs.ENSEMBL[i]
    return(colSums(gene.first.ENSEMBL==gen & essentail.single.gMCS))
  }
  
  obtain_essential_single_gMCS_simple_class <- function(i){
    gen = genes.gMCSs.ENSEMBL.12[i]
    aux1 <- gene.first.ENSEMBL==gen & essentail.single.gMCS
    
    aux2 <- apply(aux1*1, 1, sum)
    aux3 <- aux2/ncol(aux1)
    aux4 <- order(aux3, decreasing = T)
    aux3 <- matrix(aux3, ncol = 1)
    aux5 <- data.frame(gen = gen, task = "combined", gMCS = aux4,aux3[aux4])
    colnames(aux5)[-c(1:3)] <- levels(sample.class)
    aux5 <- aux5[aux5[,4]>0,]
    
    return(aux5)
  }
  
  obtain_essential_single_gMCS_multiple_class <- function(i){
    gen = genes.gMCSs.ENSEMBL.12[i]
    # gen = genes.gMCSs.ENSEMBL[i]
    # gen = "ENSG00000131471"
    # gen = "ENSG00000117143"
    aux1 <- gene.first.ENSEMBL==gen & essentail.single.gMCS
    
    if (sum(aux1)==0){
      return(NA)
    } else {
      
      aux2 <- t(rowsum(t((aux1)*1), sample.class))
      rownames(aux2) <- 1:dim(aux2)[1]
      aux3 <- aux2[order(rowSums(aux2),decreasing = T),]
      if (sum(rowSums(aux3)>0)>1){
        aux3 <- aux3[rowSums(aux3)>0,]
      } else {
        aux3 <- aux3[1:4,]
      }
      
      aux5 <- data.frame(gen = gen, task = "combined", gMCS = rownames(aux3),aux3)
      aux5 <- aux5[rowSums(aux5[,-c(1:3)])>0,]
      
      for (col in colnames(aux3)){
        aux5[,col] <- aux5[,col]/sum(sample.class==col)
      }   
      aux5 <- aux5[order(apply(aux5[,-c(1:3)],1,mean),decreasing = T),]
      
      return(aux5)
    }
  }
  
  
  nn <- length(genes.gMCSs.ENSEMBL)
  # run essential genes
  print(paste("Calculating essential genes (",nn, ") for ",ncol(gene.exp), " samples:"))
  pb <- txtProgressBar(max = nn, style = 3)
  progress <- function(i) {setTxtProgressBar(pb, i); incProgress(1/nn)}
  opts <- list(progress = progress)
  withProgress(message = 'Calculating essential genes',
               detail = '\nThis may take a while...', {
                 results <- foreach(i=1:nn, .options.snow = opts) %dopar% obtain_essential_single_genes_gMCS(i) 
               })
  
  nn <- length(genes.gMCSs.ENSEMBL.12)
  # run essential gMCS for each gen
  print(paste0("Calculating essential gMCS for each gene (",nn, ") for ",ncol(gene.exp), " samples:"))
  pb <- txtProgressBar(max = nn, style = 3)
  progress <- function(i) {setTxtProgressBar(pb, i); incProgress(1/nn)}
  opts <- list(progress = progress)
  withProgress(message = 'Calculating essential gMCS for each gene',
               detail = '\nThis may take a while...', {
                 if (length(levels(sample.class))>1){
                   list.gMCS.essential <- foreach(i=1:nn, .options.snow = opts) %dopar% obtain_essential_single_gMCS_multiple_class(i)
                 } else {
                   list.gMCS.essential <- foreach(i=1:nn, .options.snow = opts) %dopar% obtain_essential_single_gMCS_simple_class(i) 
                 }
               })
  
  
  mat.essential.gene <- do.call(rbind,results)
  rownames(mat.essential.gene) <- genes.gMCSs.ENSEMBL
  
  list.gMCS.essential <- do.call(rbind, list.gMCS.essential)
  list.gMCS.essential <- na.omit(list.gMCS.essential)
  
  
  num.essential.gene$gen <- genes.gMCSs.ENSEMBL
  num.essential.gene$num.gMCS <- apply(gMCSs.ENSEMBL.mat[,genes.gMCSs.ENSEMBL],2,sum)
  if (length(levels(sample.class))>1) {
    num.essential.gene[,levels(sample.class)] <- t(apply(mat.essential.gene, 1, function(x){rowsum((x>0)*1, sample.class)}))
  } else {
    num.essential.gene[,levels(sample.class)] <- apply(mat.essential.gene>0, 1, sum)
  }
  
  
  ratio.essential.gene <- num.essential.gene
  for (col in levels(sample.class))
  {
    ratio.essential.gene[,col] <- ratio.essential.gene[,col]/sum(sample.class==col)
  }
  
  ## Set essential genes (gMCS of only 1 gene)
  for (col in levels(sample.class))
  {
    num.essential.gene[genes.ENSEMBL.essential, col] <- sum(sample.class==col)
    ratio.essential.gene[genes.ENSEMBL.essential, col] <- 1
  }
  mat.essential.gene[genes.ENSEMBL.essential, ] <- 1
  
  
  # close cluster
  stopCluster(cl = cl)
  
  
  # Add the gene information (symbol, entrez, isEssential)
  colnames(num.essential.gene)[colnames(num.essential.gene)=="gen"] <- "ENSEMBL"
  colnames(ratio.essential.gene)[colnames(ratio.essential.gene)=="gen"] <- "ENSEMBL"
  colnames(list.gMCS.essential)[colnames(list.gMCS.essential)=="gen"] <- "ENSEMBL"
  
  num.essential.gene <- merge(table.genes.HumanGEM, num.essential.gene)
  ratio.essential.gene <- merge(table.genes.HumanGEM, ratio.essential.gene)
  list.gMCS.essential <- merge(table.genes.HumanGEM, list.gMCS.essential)
  
  
  # Export the results 
  return(list(gene.ratio = gene.ratio,
              num.essential.gene = num.essential.gene, 
              ratio.essential.gene = ratio.essential.gene,
              mat.essential.gene = mat.essential.gene,
              list.gMCS.essential = list.gMCS.essential))
  
}
