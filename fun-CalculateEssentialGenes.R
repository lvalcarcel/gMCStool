CalculateEssentialGenes <- function(gene.exp, # gene expression
                                    gMCS.info, # gMCS information
                                    sample.class, sample.cohort = NULL, #metadata of the samples
                                    threshold_logFC = 1e-3, # ratio threshold, to be used in gmcsTH y en singleTH
                                    threshold_exp = ifelse(min(gene.exp, na.rm = T)<0, -Inf, 1e-3), # minimum expression to be taken into account for the quantile calculation, only used in gmcsTH
                                    isShiny = F, # to use the tool inside and outside the shiny app
                                    parallel.nCores = 1, # if gmcsTH, decide to use parallel processing or not
                                    parallel.mode =  c("SOCK", "FORK", "FUTURE_CLUSTER", "FUTURE_MULTICORE", "FUTURE_MULTISESSION"),
                                    thresholding_methodology = c("gmcsTH", "localT2", "singleTH", 'fastcormics'), # How to calculate the threshold
                                    plot.fastcormics = F, plot.fastcormics.dir = ".", # print fastcormics information
                                    log2.fastcormics = T, # should we consider log2 in fastcormics distribution calculation?
                                    gmcsTH_perc = NULL, singleTH = NULL, # data to calculate threshold
                                    genes.in.set = rownames(gene.exp), # information to limit the genes used in localT2
                                    gene.ON.OFF = NULL, gene.ratio = NULL, localT2 = NULL, # data if threshold has been calculated and is given
                                    calculateDoubleKO = F, # Should the code calculate single or double KO?
                                    CompleteResultsDoubleKO = T, # include the list of KO with gMCS. If False, only the matrix and summary table are calculated
                                    CompleteResultsSimpleKO = T, # include the list of KO with gMCS. If False, only the matrix and summary table are calculated
                                    SimplifyGMCS = 0 # Simplify database of gMCSs
){
  
  
  # examine the inputs  
  thresholding_methodology = match.arg(thresholding_methodology)
  parallel.mode = match.arg(parallel.mode)
  if (is.null(gmcsTH_perc) & thresholding_methodology == "gmcsTH") {stop("You must provide a quantile for the threshold: gmcsTH_perc")}
  if (is.null(singleTH) & thresholding_methodology == "singleTH") {stop("You must provide a single value for the threshold: singleTH")}
  if (CompleteResultsDoubleKO) CompleteResultsSimpleKO = T
  if (!calculateDoubleKO) CompleteResultsDoubleKO = F
  if (grepl("win", Sys.info()['sysname'], ignore.case = T) & parallel.mode=="FORK") {parallel.nCores = 1} # if windows, select 1 core}
  if (length(parallel.nCores)==1){ parallel.nCores.1 = parallel.nCores.2 = parallel.nCores} else {parallel.nCores.1 = parallel.nCores[1]; parallel.nCores.2 = parallel.nCores[2]}
  if (parallel.mode!="FORK") {parallel.nCores.2 = 1} # only select one core for second part if use something different as FORK
  
  
  # Calculate first and second most expressed gene in each gMCS in each sample ####
  gene.exp[is.na(gene.exp)] <- 0 # eliminate nans
  
  
  ################################################################
  ####        s Simplification of gMCSs along tasks           ####
  ################################################################
  
  if (SimplifyGMCS>0){
    gMCS.info <- FunSimplifyGMCS(gMCS.info, trace = T, complete = SimplifyGMCS==2)
  } else {
    rownames(gMCS.info$gMCSs.ENSEMBL.mat) <- 1:nrow(gMCS.info$gMCSs.ENSEMBL.mat)
  }
  
  
  ################################################################
  ####        define highly and lowly expressed genes         ####
  ################################################################
  
  # source("fun-CalculateGenesOnOff.R")
  
  if (is.null(gene.ON.OFF)){
    # only calculate if the data is not provided
    GeneResults <- CalculateGenesOnOff(gene.exp = gene.exp, # gene expression
                                       gMCS.info = gMCS.info, # gMCS information
                                       sample.class = sample.class,
                                       sample.cohort = sample.cohort, #metadata of the samples
                                       threshold_exp = threshold_exp, # minimum expression to be taken into account for the quantile calculation, only used in gmcsTH
                                       isShiny = isShiny, # to use the tool inside and outside the shiny app
                                       parallel.nCores = parallel.nCores.1, # if gmcsTH, decide to use parallel processing or not
                                       parallel.mode =  parallel.mode,
                                       thresholding_methodology = thresholding_methodology, # How to calculate the threshold
                                       gmcsTH_perc = gmcsTH_perc,
                                       singleTH = singleTH,
                                       genes.in.set = genes.in.set) # data to calculate threshold
    
    # extract the information
    localT2 <- GeneResults$localT2
    gene.ON.OFF <- as.matrix(GeneResults$gene.ON.OFF)
    gene.ratio <- GeneResults$gene.ratio
  }
  
  
  ################################################################
  ####            extract info from the gMCSs                 ####
  ################################################################
  
  gMCSs.ENSEMBL.mat <- gMCS.info$gMCSs.ENSEMBL.mat        # matrix that relates gMCS and genes
  genes.gMCSs.ENSEMBL <- gMCS.info$genes.gMCSs.ENSEMBL    # all genes contained in the gMCS
  gMCSs.ENSEMBL.txt <- gMCS.info$gMCSs.ENSEMBL.txt        # array of unique gMCSs along tasks, text. genes separated by "--"
  if ("gMCSs.ENSEMBL.txt.list.num" %in% names(gMCS.info)){
    # gMCSs.ENSEMBL.txt.list <- gMCS.info$gMCSs.ENSEMBL.txt.list
    gMCSs.ENSEMBL.txt.list.num <- gMCS.info$gMCSs.ENSEMBL.txt.list.num
    # gMCSs.ENSEMBL.txt.list <- lapply(gMCSs.ENSEMBL.txt.list.num, function(x) genes.gMCSs.ENSEMBL[x])
  } else { 
    gMCSs.ENSEMBL.txt.list <- lapply(strsplit(gMCSs.ENSEMBL.txt, '--'),sort) # array of unique gMCSs along task, but separated
    gMCSs.ENSEMBL.txt.list.num <- lapply(gMCSs.ENSEMBL.txt.list, match, genes.gMCSs.ENSEMBL) # array of unique gMCSs along task, but separated
  }

  table.genes.HumanGEM <- gMCS.info$table.genes.HumanGEM  # table that relates genes in ENSEMBL, SYMBOL and ENTREZ_ID
  
  # ensure that it contains all the genes in the desired order for the calculation
  gene.ON.OFF <- gene.ON.OFF[genes.gMCSs.ENSEMBL,]
  
  ################################################################
  ####        single and double essential genes               ####
  ################################################################
  
  # browser()
  index_gmcs_length_1 <- as.numeric(which(rowSums(gMCSs.ENSEMBL.mat)==1))
  
  genes.ENSEMBL.essential <- gMCSs.ENSEMBL.txt[index_gmcs_length_1]
  genes.ENSEMBL.essential.num <- match(genes.ENSEMBL.essential, genes.gMCSs.ENSEMBL)
  
  if (calculateDoubleKO){
    index_gmcs_length_2 <- as.numeric(which(rowSums(gMCSs.ENSEMBL.mat)==2))
    # genes.ENSEMBL.essential.pair <- lapply(gMCSs.ENSEMBL.txt.list[index_gmcs_length_2], sort)
    genes.ENSEMBL.essential.pair <- lapply(gMCSs.ENSEMBL.txt.list.num[index_gmcs_length_2], function(x) sort(genes.gMCSs.ENSEMBL[x]))
  }
  
  # define a set for single gene KO ### 
  
  num.essential.gene <- as.data.frame(matrix(0,nrow = length(genes.gMCSs.ENSEMBL),ncol = length(levels(sample.class))+2))
  colnames(num.essential.gene) <- c("gen","num.gMCS",levels(sample.class))
  rownames(num.essential.gene) <- genes.gMCSs.ENSEMBL
  num.essential.gene$gen <- genes.gMCSs.ENSEMBL
  
  mat.essential.gene <- Matrix(0,nrow = length(genes.gMCSs.ENSEMBL),ncol = dim(gene.exp)[2])
  colnames(mat.essential.gene) <- colnames(gene.exp)
  rownames(mat.essential.gene) <- 1:length(genes.gMCSs.ENSEMBL)
  
  if (CompleteResultsSimpleKO & !isShiny){
    list.gene.of.gMCS.by.sample <- Matrix(0,nrow = nrow(gMCSs.ENSEMBL.mat),ncol = dim(gene.exp)[2])
    list.gene.of.gMCS.by.sample <- lapply(genes.gMCSs.ENSEMBL, function(x){ return(list.gene.of.gMCS.by.sample)})
    names(list.gene.of.gMCS.by.sample) <- genes.gMCSs.ENSEMBL
  } else {
    list.gene.of.gMCS.by.sample <- F
  }
  list.gene.of.gMCS.by.sample.aux <- list() # auxiliary list
  
  if (calculateDoubleKO){ 
    
    # browser()
    # define a set for double gene KO ###
    
    # genes.double.gMCSs.ENSEMBL <- lapply(gMCSs.ENSEMBL.txt.list[-index_gmcs_length_1], sort)
    genes.double.gMCSs.ENSEMBL <- lapply(gMCSs.ENSEMBL.txt.list.num[-index_gmcs_length_1], function(x) sort(genes.gMCSs.ENSEMBL[x]))
    
    # genes.double.gMCSs.ENSEMBL <- lapply(genes.double.gMCSs.ENSEMBL, sort)
    genes.double.gMCSs.ENSEMBL <- lapply(genes.double.gMCSs.ENSEMBL, gRbase::combn_prim, 2)
    genes.double.gMCSs.ENSEMBL <- t(do.call(cbind, genes.double.gMCSs.ENSEMBL))
    
    genes.double.gMCSs.ENSEMBL.counts <- data.frame(gen1 = genes.double.gMCSs.ENSEMBL[,1], gen2 = genes.double.gMCSs.ENSEMBL[,2]) %>% count(gen1, gen2)
    
    MagicNumber <- 10^ceiling(log10(length(genes.gMCSs.ENSEMBL))) # this magic number allow us to transform the gene pair into a single number, being g1_idx*MagicNumber + g2_idx
    
    genes.double.gMCSs.ENSEMBL <- genes.double.gMCSs.ENSEMBL.counts[,1:2]
    genes.double.gMCSs.ENSEMBL[,3] <- match(genes.double.gMCSs.ENSEMBL[,1], genes.gMCSs.ENSEMBL)
    genes.double.gMCSs.ENSEMBL[,4] <- match(genes.double.gMCSs.ENSEMBL[,2], genes.gMCSs.ENSEMBL)
    genes.double.gMCSs.ENSEMBL[,5] <- apply(genes.double.gMCSs.ENSEMBL[,1:2], 1, paste, collapse = "_")
    genes.double.gMCSs.ENSEMBL[,6] <- genes.double.gMCSs.ENSEMBL[,3]*MagicNumber + genes.double.gMCSs.ENSEMBL[,4]
    genes.double.gMCSs.ENSEMBL[,7] <- genes.double.gMCSs.ENSEMBL[,4]*MagicNumber + genes.double.gMCSs.ENSEMBL[,3]
    

    num.essential.pair.gene <- as.data.frame(matrix(0,nrow = nrow(genes.double.gMCSs.ENSEMBL),ncol = length(levels(sample.class))+5))
    colnames(num.essential.pair.gene) <- c("gen1","gen2", "essentialGene", "essentialPair", "num.gMCS",levels(sample.class))
    num.essential.pair.gene$gen1 <- genes.double.gMCSs.ENSEMBL.counts[,1]
    num.essential.pair.gene$gen2 <- genes.double.gMCSs.ENSEMBL.counts[,2]
    num.essential.pair.gene$num.gMCS <- genes.double.gMCSs.ENSEMBL.counts[,3]
    
    mat.essential.pair.gene <- Matrix(0,nrow = nrow(genes.double.gMCSs.ENSEMBL),ncol = dim(gene.exp)[2])
    colnames(mat.essential.pair.gene) <- colnames(gene.exp)
    rownames(mat.essential.pair.gene) <- 1:nrow(mat.essential.pair.gene)
    
    if (CompleteResultsDoubleKO & !isShiny){
      list.pairgene.of.gMCS.by.sample <- Matrix(0,nrow = nrow(gMCSs.ENSEMBL.mat),ncol = dim(gene.exp)[2])
      list.pairgene.of.gMCS.by.sample <- lapply(rownames(mat.essential.pair.gene) , function(x){ return(list.pairgene.of.gMCS.by.sample)})
      names(list.pairgene.of.gMCS.by.sample) <- rownames(mat.essential.pair.gene)
    } else {
      list.pairgene.of.gMCS.by.sample <- NULL
    }
    list.pairgene.of.gMCS.by.sample.aux <- list()  # auxiliary list
  }
  
  
  # mat.essential.gene.2 <- mat.essential.gene
  # mat.essential.pair.gene.2 <- mat.essential.pair.gene
  
  # browser()
  genes.in.gmcs.over.th <- gMCSs.ENSEMBL.mat[,genes.gMCSs.ENSEMBL] %*% as.matrix(gene.ON.OFF)*1

  # store as a binary matrix, and transpose (save time for calculations)
  gMCSs.ENSEMBL.mat.binary <- t(gMCSs.ENSEMBL.mat>0.5)
  
  ## Generate functions to calculate essential genes
  fun_single_KO <- function(i){
    # gmcs with only one active gene
    idx.1 <- which(genes.in.gmcs.over.th[,i]==1, useNames = FALSE)  
    # genes that are highly expressed for the patient
    gene.up.patient <- which(gene.ON.OFF[,i], useNames = FALSE)  
    # the intersection inform us of the active genes in the gmcs with highly expressed genes
    idx.3 <- which(gMCSs.ENSEMBL.mat.binary[gene.up.patient, idx.1], arr.ind = T, useNames = FALSE)
    idx.3 <- gene.up.patient[idx.3[,1]]
    # essentiality matrixes
    return(data.frame(gMCS = idx.1,
                      ENSEMBL = idx.3,
                      sample = i, 
                      stringsAsFactors = F))
  }
  
  
  fun_double_KO <- function(i){
    # gmcs with only two active genes
    idx.2 <- which(genes.in.gmcs.over.th[,i]==2, useNames = FALSE)  
    # genes that are highly expressed for the patient
    gene.up.patient <- which(gene.ON.OFF[,i], useNames = FALSE)  
    # the intersection inform us of the active genes in the gmcs with highly expressed genes
    idx.4 <- which(gMCSs.ENSEMBL.mat.binary[gene.up.patient, idx.2], arr.ind = T, useNames = FALSE)
    idx.4 <- matrix(idx.4[,1], nrow = 2, byrow = F)
    idx.4 <- matrix(c(gene.up.patient[idx.4[1,]], gene.up.patient[idx.4[2,]]), ncol = 2, byrow = F)
    # remove genes that are considered as essential for that sample
    idx.aux <- ((idx.4[,1] %in% list.gene.of.gMCS.by.sample.aux[[i]]$ENSEMBL) | (idx.4[,2] %in% list.gene.of.gMCS.by.sample.aux[[i]]$ENSEMBL))
    idx.aux <- idx.aux | ((idx.4[,1] %in% genes.ENSEMBL.essential.num) | (idx.4[,2] %in% genes.ENSEMBL.essential.num))
    idx.aux <- !(idx.aux)
    # Transform into double gene index
    # idx.4 <- apply(idx.4, 1, paste, collapse = "_")
    idx.4 <- idx.4[,1]*MagicNumber + idx.4[,2]
    idx.4 <- matrix(c(match(idx.4, genes.double.gMCSs.ENSEMBL[,6]), 
                      match(idx.4, genes.double.gMCSs.ENSEMBL[,7])),
                    ncol = 2, byrow = F, nrow = length(idx.4))
    idx.4 <- rowSums(idx.4, na.rm = T) # faster than apply(x,1,sum,na.rm = T)
    # essentiality matrixes
    # here we store the gmcs that are responsible of the essentiality
    return(data.frame(gMCS = idx.2[idx.aux],
                      ENSEMBL = idx.4[idx.aux],
                      sample = i, 
                      stringsAsFactors = F))
  }
  
  
  
  # browser()
  # function to search for essential genes (multicore or singlecore)
  if (parallel.nCores.2<=1){
    print(paste("Calculating essential gMCS for each gene  for",ncol(gene.exp), "samples:"))
    
    pb <- txtProgressBar(max = ncol(gene.ON.OFF), style = 3)
    if (isShiny){
      withProgress(message = 'Calculating essential genes',
                   detail = '\nThis may take a while...', {
                     for (i in 1:ncol(gene.ON.OFF)){
                       setTxtProgressBar(pb, i); incProgress(1/ncol(gene.ON.OFF))
                       list.gene.of.gMCS.by.sample.aux[[i]] <- fun_single_KO(i)
                     }                   })
    } else {
      for (i in 1:ncol(gene.ON.OFF)){
        setTxtProgressBar(pb, i)
        list.gene.of.gMCS.by.sample.aux[[i]] <- fun_single_KO(i)
      }
    }
    close(pb)
    
    # repeat for doble KO
    if (calculateDoubleKO){
      print(paste("Calculating essential gMCS for each pair of genes for",ncol(gene.exp), "samples:"))
      
      pb <- txtProgressBar(max = ncol(gene.ON.OFF), style = 3)
      if (isShiny){
        withProgress(message = 'Calculating essential gene pairs',
                     detail = '\nThis may take a while...', {
                       for (i in 1:ncol(gene.ON.OFF)){
                         setTxtProgressBar(pb, i); incProgress(1/ncol(gene.ON.OFF))
                         list.pairgene.of.gMCS.by.sample.aux[[i]] <- fun_double_KO(i)
                       }
                     })
      } else {
        for (i in 1:ncol(gene.ON.OFF)){
          setTxtProgressBar(pb, i)
          list.pairgene.of.gMCS.by.sample.aux[[i]] <- fun_double_KO(i)
        }
      }
      close(pb)
    }
  } else if (parallel.mode == "FORK"){
    
    print(paste("Calculating essential gMCS using multicore processing with", parallel.nCores.2, "cores for each gene for",ncol(gene.exp), "samples:"))
    list.gene.of.gMCS.by.sample.aux <- pbmcapply::pbmclapply(1:ncol(gene.ON.OFF),
                                                             fun_single_KO, 
                                                             mc.cores = parallel.nCores.2,
                                                             mc.preschedule = T,
                                                             mc.style = "txt", mc.substyle = 3)
    if (calculateDoubleKO){
      
      print(paste("Calculating essential gMCS using multicore processing with", parallel.nCores.2, "cores for each pair of genes for",ncol(gene.exp), "samples:"))
      list.pairgene.of.gMCS.by.sample.aux <- pbmcapply::pbmclapply(1:ncol(gene.ON.OFF),
                                                                   fun_double_KO,
                                                                   mc.cores = parallel.nCores.2,
                                                                   mc.preschedule = T,
                                                                   mc.style = "txt", mc.substyle = 3)
    }
    
    
  } 
  
  print("Processing the results....")
  pb <- txtProgressBar(max = ifelse(calculateDoubleKO, 13, 7), style = 3); pb_cont <- 1
  
  # browser()
  # Fill the matrix
  list.gene.of.gMCS.by.sample.aux <- do.call(rbind, list.gene.of.gMCS.by.sample.aux)
  
  aux <- list.gene.of.gMCS.by.sample.aux %>% dplyr::select(ENSEMBL, sample) %>% unique() %>% as.matrix()
  mat.essential.gene[aux] <- 1
  
  list.gene.of.gMCS.by.sample.aux <- split(list.gene.of.gMCS.by.sample.aux, #%>% mutate(ENSEMBL = rownames(mat.essential.gene)[ENSEMBL]),
                                           ~ENSEMBL)
  setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1

  if (calculateDoubleKO){ 
    list.pairgene.of.gMCS.by.sample.aux <- do.call(rbind, list.pairgene.of.gMCS.by.sample.aux)
    
    aux <- list.pairgene.of.gMCS.by.sample.aux %>% dplyr::select(ENSEMBL, sample) %>% unique() %>% as.matrix()
    mat.essential.pair.gene[aux] <- 1
    
    list.pairgene.of.gMCS.by.sample.aux <- split(list.pairgene.of.gMCS.by.sample.aux, #%>% mutate(ENSEMBL = rownames(mat.essential.gene)[ENSEMBL]),
                                                 ~ENSEMBL)

    setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
  }
  
  if (CompleteResultsSimpleKO & !isShiny){
    for (g in names(list.gene.of.gMCS.by.sample.aux)){
      aux <- list.gene.of.gMCS.by.sample.aux[[g]]
      list.gene.of.gMCS.by.sample[[g]][cbind(aux$gMCS,aux$sample)] <- rep(1, nrow(aux))
    }
  }
  setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
  
  
  if (CompleteResultsDoubleKO & !isShiny){
    list.pairgene.of.gMCS.by.sample.aux <- do.call(rbind, list.pairgene.of.gMCS.by.sample.aux)
    list.pairgene.of.gMCS.by.sample.aux <- split(list.pairgene.of.gMCS.by.sample.aux, # %>%  mutate(ENSEMBL = genes.double.gMCSs.ENSEMBL[ENSEMBL,5]),
                                                 ~ENSEMBL)
    
    for (g in names(list.pairgene.of.gMCS.by.sample.aux)){
      aux <- list.pairgene.of.gMCS.by.sample.aux[[g]]
      list.pairgene.of.gMCS.by.sample[[g]][cbind(aux$gMCS,aux$sample)] <- rep(1, nrow(aux))
    }
  }
  setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
  
  if (CompleteResultsSimpleKO){
    
    list.gMCS.essential <- reshape2::dcast(gen + gMCS ~ class,
                                           data = do.call(rbind, list.gene.of.gMCS.by.sample.aux) %>%
                                             mutate(gen = ENSEMBL) %>%
                                             mutate(class = sample.class[sample]),
                                           fun.aggregate = length, value.var = "class")
    list.gMCS.essential[,levels(sample.class)] <- t(apply(list.gMCS.essential[,levels(sample.class)],1,function(x){x/table(sample.class)}, simplify = T))
    list.gMCS.essential <- list.gMCS.essential %>% mutate(gen = genes.gMCSs.ENSEMBL[as.numeric(as.character(gen))])
    
  } else {
    list.gMCS.essential <- NULL
  }
  setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
  
  # generate summary tables
  rownames(mat.essential.gene) <- genes.gMCSs.ENSEMBL
  num.essential.gene$gen <- genes.gMCSs.ENSEMBL
  num.essential.gene$num.gMCS <- colSums(gMCSs.ENSEMBL.mat[,genes.gMCSs.ENSEMBL])
  
  if (length(levels(sample.class))>1) {
    num.essential.gene[,levels(sample.class)] <- t(rowsum(t(as.matrix(mat.essential.gene>0)*1), sample.class))
  } else {
    num.essential.gene[,levels(sample.class)] <- rowSums(mat.essential.gene>0)
  }
  setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
  
  ratio.essential.gene <- num.essential.gene
  for (col in levels(sample.class)){
    ratio.essential.gene[,col] <- ratio.essential.gene[,col]/sum(sample.class==col)
  }
  
  ## Set essential genes (gMCS of only 1 gene)
  for (col in levels(sample.class)){
    num.essential.gene[genes.ENSEMBL.essential, col] <- sum(sample.class==col)
    ratio.essential.gene[genes.ENSEMBL.essential, col] <- 1
  }
  mat.essential.gene[genes.ENSEMBL.essential, ] <- 1
  setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
  
  # Add the gene information (symbol, entrez, isEssential)
  colnames(num.essential.gene)[colnames(num.essential.gene)=="gen"] <- "ENSEMBL"
  colnames(ratio.essential.gene)[colnames(ratio.essential.gene)=="gen"] <- "ENSEMBL"
  
  num.essential.gene <- merge(table.genes.HumanGEM, num.essential.gene)
  ratio.essential.gene <- merge(table.genes.HumanGEM, ratio.essential.gene)
  
  if (CompleteResultsSimpleKO){
    colnames(list.gMCS.essential)[colnames(list.gMCS.essential)=="gen"] <- "ENSEMBL"
    list.gMCS.essential <- merge(table.genes.HumanGEM, list.gMCS.essential)
  }
  setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
  
  # analyze results if doble KO is calculated
  if (calculateDoubleKO){
    
    # matrix of essentiality due to the double KO
    mat.essential.single.pair.combined.gene <- mat.essential.pair.gene
    
    # matrix should include the single KOs
    mat.essential.single.pair.combined.gene <- mat.essential.single.pair.combined.gene + mat.essential.gene[genes.double.gMCSs.ENSEMBL[,3],]
    mat.essential.single.pair.combined.gene <- mat.essential.single.pair.combined.gene + mat.essential.gene[genes.double.gMCSs.ENSEMBL[,4],]
    
    if (CompleteResultsDoubleKO) {
      list.gMCS.essential.pair <- reshape2::dcast(genes + gMCS ~ class,
                                                  data = do.call(rbind, list.pairgene.of.gMCS.by.sample.aux) %>% 
                                                    mutate(class = sample.class[sample]) %>% 
                                                    mutate(genes = ENSEMBL),
                                                  fun.aggregate = length, value.var = "class")
      list.gMCS.essential.pair[,levels(sample.class)] <- t(apply(list.gMCS.essential.pair[,levels(sample.class)],1,function(x){x/table(sample.class)}, simplify = T))
      list.gMCS.essential.pair <- list.gMCS.essential.pair %>% 
        mutate(genes = genes.double.gMCSs.ENSEMBL[as.numeric(as.character(genes)),5]) %>% 
        mutate(gen1 = sub("_.*", "",genes)) %>% 
        mutate(gen2 = sub(".*_", "",genes)) %>% 
        mutate(gMCS = as.character(gMCS))
      list.gMCS.essential.pair$gMCS <- rownames(gMCSs.ENSEMBL.mat)[as.numeric(as.character(list.gMCS.essential.pair$gMCS))]
      list.gMCS.essential.pair <- list.gMCS.essential.pair[,c("gen1", "gen2", "gMCS", levels(sample.class))]
      
      
    } else {list.gMCS.essential.pair <- NULL}
    setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
    
    # add gene pair name
    rownames(mat.essential.pair.gene) <- genes.double.gMCSs.ENSEMBL[,5]
    rownames(mat.essential.single.pair.combined.gene) <- genes.double.gMCSs.ENSEMBL[,5]
    
    # generate summary tables
    num.essential.pair.gene$essentialPair <- apply(num.essential.pair.gene[,c("gen1", "gen2")], 1, paste0, collapse = "_") %in%
      unlist(lapply(genes.ENSEMBL.essential.pair, paste0, collapse = "_"))
    
    num.essential.pair.gene$essentialGene <- (num.essential.pair.gene$gen1 %in% genes.ENSEMBL.essential)*1 + (num.essential.pair.gene$gen2 %in% genes.ENSEMBL.essential)*2
    
    num.essential.single.pair.combined.gene <- num.essential.pair.gene
    
    if (length(levels(sample.class))>1) {
      num.essential.pair.gene[,levels(sample.class)] <- t(rowsum(t(as.matrix(mat.essential.pair.gene>0)*1), sample.class))
      num.essential.single.pair.combined.gene[,levels(sample.class)] <- t(rowsum(t(as.matrix(mat.essential.single.pair.combined.gene>0)*1), sample.class))
    } else {
      num.essential.pair.gene[,levels(sample.class)] <- rowSums(mat.essential.pair.gene>0)
      num.essential.single.pair.combined.gene[,levels(sample.class)] <- rowSums(mat.essential.single.pair.combined.gene>0)
    }
    setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
    
    
    ratio.essential.pair.gene <- num.essential.pair.gene
    ratio.essential.single.pair.combined.gene <- num.essential.single.pair.combined.gene
    for (col in levels(sample.class)){
      ratio.essential.pair.gene[,col] <- num.essential.pair.gene[,col]/sum(sample.class==col)
      ratio.essential.single.pair.combined.gene[,col] <- num.essential.single.pair.combined.gene[,col]/sum(sample.class==col)
    }
    
    
    ## Set essential genes (gMCS of only 2 gene)
    for (col in levels(sample.class)){
      num.essential.pair.gene[sapply(genes.ENSEMBL.essential.pair,paste, collapse = "_"), col] <- sum(sample.class==col)
      ratio.essential.pair.gene[sapply(genes.ENSEMBL.essential.pair,paste, collapse = "_"), col] <- 1
    }
    mat.essential.pair.gene[sapply(genes.ENSEMBL.essential.pair,paste, collapse = "_"), ] <- 1
    setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
    
    # Add the gene information (symbol, isEssentialPair)
    
    
    colnames(num.essential.pair.gene)[colnames(num.essential.pair.gene)=="gen1"] <- "ENSEMBL.1"
    colnames(num.essential.pair.gene)[colnames(num.essential.pair.gene)=="gen2"] <- "ENSEMBL.2"
    
    # num.essential.pair.gene.2 <- num.essential.pair.gene
    num.essential.pair.gene <- merge(table.genes.HumanGEM %>% as.data.frame() %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.2 = ENSEMBL) %>% dplyr::rename(SYMBOL.2 = SYMBOL),
                                     num.essential.pair.gene %>% as.data.frame())
    num.essential.pair.gene <- merge(table.genes.HumanGEM %>% as.data.frame() %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.1 = ENSEMBL) %>% dplyr::rename(SYMBOL.1 = SYMBOL),
                                     num.essential.pair.gene %>% as.data.frame())
    
    colnames(ratio.essential.pair.gene)[colnames(ratio.essential.pair.gene)=="gen1"] <- "ENSEMBL.1"
    colnames(ratio.essential.pair.gene)[colnames(ratio.essential.pair.gene)=="gen2"] <- "ENSEMBL.2"
    
    ratio.essential.pair.gene <- merge(table.genes.HumanGEM %>%  dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.2 = ENSEMBL) %>% dplyr::rename(SYMBOL.2 = SYMBOL),
                                       ratio.essential.pair.gene)
    ratio.essential.pair.gene <- merge(table.genes.HumanGEM %>%  dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.1 = ENSEMBL) %>% dplyr::rename(SYMBOL.1 = SYMBOL),
                                       ratio.essential.pair.gene)
    
    if (CompleteResultsDoubleKO){
      colnames(list.gMCS.essential.pair)[colnames(list.gMCS.essential.pair)=="gen1"] <- "ENSEMBL.1"
      colnames(list.gMCS.essential.pair)[colnames(list.gMCS.essential.pair)=="gen2"] <- "ENSEMBL.2"
      
      list.gMCS.essential.pair <- merge(table.genes.HumanGEM %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.2 = ENSEMBL) %>% dplyr::rename(SYMBOL.2 = SYMBOL),
                                        list.gMCS.essential.pair)
      list.gMCS.essential.pair <- merge(table.genes.HumanGEM %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.1 = ENSEMBL) %>% dplyr::rename(SYMBOL.1 = SYMBOL),
                                        list.gMCS.essential.pair)
    }
    setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
    
    colnames(num.essential.single.pair.combined.gene)[colnames(num.essential.single.pair.combined.gene)=="gen1"] <- "ENSEMBL.1"
    colnames(num.essential.single.pair.combined.gene)[colnames(num.essential.single.pair.combined.gene)=="gen2"] <- "ENSEMBL.2"
    
    num.essential.single.pair.combined.gene <- merge(table.genes.HumanGEM %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.2 = ENSEMBL) %>% dplyr::rename(SYMBOL.2 = SYMBOL),
                                                     num.essential.single.pair.combined.gene)
    num.essential.single.pair.combined.gene <- merge(table.genes.HumanGEM %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.1 = ENSEMBL) %>% dplyr::rename(SYMBOL.1 = SYMBOL),
                                                     num.essential.single.pair.combined.gene)
    
    colnames(ratio.essential.single.pair.combined.gene)[colnames(ratio.essential.single.pair.combined.gene)=="gen1"] <- "ENSEMBL.1"
    colnames(ratio.essential.single.pair.combined.gene)[colnames(ratio.essential.single.pair.combined.gene)=="gen2"] <- "ENSEMBL.2"
    
    ratio.essential.single.pair.combined.gene <- merge(table.genes.HumanGEM %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.2 = ENSEMBL) %>% dplyr::rename(SYMBOL.2 = SYMBOL),
                                                       ratio.essential.single.pair.combined.gene)
    ratio.essential.single.pair.combined.gene <- merge(table.genes.HumanGEM %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::rename(ENSEMBL.1 = ENSEMBL) %>% dplyr::rename(SYMBOL.1 = SYMBOL),
                                                       ratio.essential.single.pair.combined.gene)
    setTxtProgressBar(pb, pb_cont); pb_cont <- pb_cont + 1
    
  } else {
    num.essential.pair.gene = NULL
    ratio.essential.pair.gene = NULL
    num.essential.single.pair.combined.gene = NULL
    ratio.essential.single.pair.combined.gene = NULL
    mat.essential.pair.gene = NULL
    list.gMCS.essential.pair = NULL
    list.pairgene.of.gMCS.by.sample = NULL
  }
  close(pb)
  
  # print(pb_cont)
  
  # Export the results
  return(list(gene.ratio = gene.ratio,
              gene.ON.OFF = gene.ON.OFF,
              localT2 = localT2,
              num.essential.gene = num.essential.gene,
              ratio.essential.gene = ratio.essential.gene,
              mat.essential.gene = mat.essential.gene,
              list.gene.of.gMCS.by.sample = list.gene.of.gMCS.by.sample,
              list.gMCS.essential = list.gMCS.essential,
              num.essential.pair.gene = num.essential.pair.gene,
              ratio.essential.pair.gene = ratio.essential.pair.gene,
              num.essential.single.pair.combined.gene = num.essential.single.pair.combined.gene,
              ratio.essential.single.pair.combined.gene = ratio.essential.single.pair.combined.gene,
              mat.essential.pair.gene = mat.essential.pair.gene,
              list.gMCS.essential.pair = list.gMCS.essential.pair,
              list.pairgene.of.gMCS.by.sample = list.pairgene.of.gMCS.by.sample ))
  
  
}
