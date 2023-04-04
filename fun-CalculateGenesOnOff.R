CalculateGenesOnOff <- function(gene.exp, # gene expression
                                gMCS.info, # gMCS information
                                sample.class, sample.cohort = NULL, #metadata of the samples
                                isShiny = F, # to use the tool inside and outside the shiny app
                                parallel.nCores = 1, # if gmcsTH, decide to use parallel processing or not
                                parallel.mode =  c("SOCK", "FORK", "FUTURE_CLUSTER", "FUTURE_MULTICORE", "FUTURE_MULTISESSION"),
                                thresholding_methodology = c("gmcsTH", "localT2", "singleTH", 'fastcormics'), # How to calculate the threshold
                                gmcsTH_perc = NULL, singleTH = NULL, # data to calculate threshold
                                plot.fastcormics = F, plot.fastcormics.dir = ".", # print fastcormics information
                                log2.fastcormics = T, # should we consider log2 in fastcormics distribution calculation?
                                genes.in.set = rownames(gene.exp), # information to limit the genes used in localT2
                                threshold_logFC = 1e-3, # ratio threshold, to be used in gmcsTH and singleTH
                                threshold_exp = ifelse(min(gene.exp)<0, -Inf, 1e-3) # minimum expression to be taken into account for the quantile calculation, only used in gmcsTH
){
  
  
  # examine the inputs
  thresholding_methodology = match.arg(thresholding_methodology)
  parallel.mode = match.arg(parallel.mode)
  if (is.null(gmcsTH_perc) & thresholding_methodology == "gmcsTH") {stop("You must provide a quantile for the threshold: gmcsTH_perc")}
  if (is.null(singleTH) & thresholding_methodology == "singleTH") {stop("You must provide a single value for the threshold: singleTH")}
  if (grepl("win", Sys.info()['sysname'], ignore.case = T) & parallel.mode=="FORK") {parallel.nCores = 1} # if windows, select 1 core}
  # if (is.null(threshold_exp)) {threshold_exp <- ifelse(min(gene.exp)<0, -Inf, 1e-3)}
  
  # extract info from the gMCSs
  gMCSs.ENSEMBL.mat <- gMCS.info$gMCSs.ENSEMBL.mat        # matrix that relates gMCS and genes
  genes.gMCSs.ENSEMBL <- gMCS.info$genes.gMCSs.ENSEMBL    # all genes contained in the gMCS
  gMCSs.ENSEMBL.txt <- gMCS.info$gMCSs.ENSEMBL.txt        # array of unique gMCSs along tasks, text. genes separated by "--"
  if ("gMCSs.ENSEMBL.txt.list.num" %in% names(gMCS.info)){
    # gMCSs.ENSEMBL.txt.list <- gMCS.info$gMCSs.ENSEMBL.txt.list
    gMCSs.ENSEMBL.txt.list.num <- gMCS.info$gMCSs.ENSEMBL.txt.list.num
  } else if (thresholding_methodology == "gmcsTH"){
    gMCSs.ENSEMBL.txt.list <- lapply(strsplit(gMCSs.ENSEMBL.txt, '--'), sort) # array of unique gMCSs along task, but separated
    gMCSs.ENSEMBL.txt.list.num <- lapply(gMCSs.ENSEMBL.txt.list, match, genes.gMCSs.ENSEMBL) # array of unique gMCSs along task, but separated
  }
  table.genes.HumanGEM <- gMCS.info$table.genes.HumanGEM  # table that relates genes in ENSEMBL, SYMBOL and ENTREZ_ID
  
  # clear names of variables to make them faster
  # names(gMCSs.ENSEMBL.txt) <- names(gMCSs.ENSEMBL.txt.list) <- names(gMCSs.ENSEMBL.txt.list.num) <- NULL
  names(gMCSs.ENSEMBL.txt) <- names(gMCSs.ENSEMBL.txt.list.num) <- NULL
  
  # Calculate first and second most expressed gene in each gMCS in each sample ###
  gene.exp[is.na(gene.exp)] <- 0 # eliminate nans
  
  ## Define the different cases
  if (thresholding_methodology == "gmcsTH"){
    ###########################  gmcsTH methodology   ###########################
    
    nn <- dim(gMCSs.ENSEMBL.mat)[1]
    
    gene.first.exp = gene.first.ENSEMBL =
      matrix(NA,nrow = nn, ncol = dim(gene.exp)[2])
    
    colnames(gene.first.exp) = colnames(gene.first.ENSEMBL) =
      colnames(gene.exp)
    
    gene.exp.gMCSs <- gene.exp[genes.gMCSs.ENSEMBL,]
    rownames(gene.exp.gMCSs) <- genes.gMCSs.ENSEMBL
    
    print(paste("Processing the",nn, "gMCSs for",ncol(gene.exp), "samples:"))
    
    
    # function to search the maximum expressed gene
    obtain_max_gene_gMCS <- function(gmcs, gMCSs.ENSEMBL.exp){
      aux <- gMCSs.ENSEMBL.exp[gmcs]
      aux <- gmcs[which.max(aux)]
      return(aux)
    }
    obtain_max_gene_gMCS_sample <- function(gMCSs.ENSEMBL.exp, gmcs.list, genes.gMCSs.ENSEMBL, obtain_max_gene_gMCS){
      print(head(gMCSs.ENSEMBL.exp))
      x <- vapply(gmcs.list, FUN = obtain_max_gene_gMCS, FUN.VALUE = numeric(1), gMCSs.ENSEMBL.exp = gMCSs.ENSEMBL.exp)
      return(matrix(genes.gMCSs.ENSEMBL[x],ncol = 1))
    }
    
    # browser()
    if (parallel.nCores<=1){
      if (isShiny){
        pb <- txtProgressBar(max = ncol(gene.exp), style = 3)
        withProgress(message = 'Calculating most expressed gene by gMCS and sample',
                     detail = '\nThis may take a while...', {
                       for (SAMPLE in 1:ncol(gene.exp)){
                         setTxtProgressBar(pb, SAMPLE); incProgress(1/ncol(gene.exp))
                         x <- gene.exp.gMCSs[,SAMPLE]
                         x <- vapply(gMCSs.ENSEMBL.txt.list.num, FUN = obtain_max_gene_gMCS, FUN.VALUE = numeric(1), gMCSs.ENSEMBL.exp = x)
                         gene.first.ENSEMBL[,SAMPLE] <- genes.gMCSs.ENSEMBL[x]
                       }
                     })
        close(pb)
      } else {
        pb <- txtProgressBar(max = ncol(gene.exp), style = 3)
        for (SAMPLE in 1:ncol(gene.exp)){
          setTxtProgressBar(pb, SAMPLE)
          x <- gene.exp.gMCSs[,SAMPLE]
          x <- vapply(gMCSs.ENSEMBL.txt.list.num, FUN = obtain_max_gene_gMCS, FUN.VALUE = numeric(1), gMCSs.ENSEMBL.exp = x)
          gene.first.ENSEMBL[,SAMPLE] <- genes.gMCSs.ENSEMBL[x]
        }
        close(pb)
      }
      
    } else if (parallel.mode == "FORK"){
      print(paste("Using multicore processing with", parallel.nCores, "cores for",ncol(gene.exp), "samples:"))
      obtain_max_gene_gMCS_sample <- function(SAMPLE){
        gMCSs.ENSEMBL.exp <- gene.exp.gMCSs[,SAMPLE]
        x <- vapply(gMCSs.ENSEMBL.txt.list.num,
                    FUN = obtain_max_gene_gMCS,
                    FUN.VALUE = numeric(1),
                    gMCSs.ENSEMBL.exp = gMCSs.ENSEMBL.exp)
        return(matrix(genes.gMCSs.ENSEMBL[x],ncol = 1))
      }
      gene.first.ENSEMBL <- pbmcapply::pbmclapply(1:ncol(gene.exp),
                                                  obtain_max_gene_gMCS_sample,
                                                  mc.cores = parallel.nCores,
                                                  mc.preschedule = T,
                                                  mc.style = "txt", mc.substyle = 3)
      
      gene.first.ENSEMBL <- do.call(cbind, gene.first.ENSEMBL)
    } else if (parallel.mode == "SOCK"){
      
      clF <- makePSOCKcluster(as.integer(parallel.nCores))
      
      obtain_max_gene_gMCS_sample <- function(gMCSs.ENSEMBL.exp){
        x <- vapply(gMCSs.ENSEMBL.txt.list.num, FUN = obtain_max_gene_gMCS, FUN.VALUE = numeric(1), gMCSs.ENSEMBL.exp = gMCSs.ENSEMBL.exp)
        return(matrix(genes.gMCSs.ENSEMBL[x],ncol = 1))
      }
      clusterExport(cl=clF, c('gMCSs.ENSEMBL.txt.list.num', 'genes.gMCSs.ENSEMBL', "obtain_max_gene_gMCS", 'obtain_max_gene_gMCS_sample'), envir = environment())
      environment(obtain_max_gene_gMCS_sample) <- .GlobalEnv
      gene.first.ENSEMBL <- ParallelLogger::clusterApply(cl = clF,
                                                         x = as.list(as.data.frame(gene.exp.gMCSs)),
                                                         fun = obtain_max_gene_gMCS_sample,
                                                         progressBar = TRUE)
      stopCluster(cl = clF)
      gene.first.ENSEMBL <- do.call(cbind, gene.first.ENSEMBL)
      
    } else if (grepl("FUTURE", parallel.mode, ignore.case = T)){
      
      if (parallel.mode=="FUTURE_CLUSTER"){
        oplan <- plan(cluster, workers = parallel.nCores)
      } else if (parallel.mode=="FUTURE_MULTICORE"){
        oplan <- plan(multicore, workers = parallel.nCores)
      } else if (parallel.mode=="FUTURE_MULTISESSION"){
        oplan <- plan(multisession  , workers = parallel.nCores)
      }
      
      obtain_max_gene_gMCS_sample <- function(gMCSs.ENSEMBL.exp, gmcs.list, genes.gMCSs.ENSEMBL, obtain_max_gene_gMCS){
        x <- vapply(gmcs.list, FUN = obtain_max_gene_gMCS, FUN.VALUE = numeric(1), gMCSs.ENSEMBL.exp = gMCSs.ENSEMBL.exp)
        return(matrix(genes.gMCSs.ENSEMBL[x],ncol = 1))
      }
      
      
      print(paste("Using multicore processing with", parallel.nCores, "cores for",ncol(gene.exp), "samples:"))
      # pb <- txtProgressBar(max = ncol(gene.exp), style = 3)
      gene.first.ENSEMBL <- future.apply::future_lapply(1:ncol(gene.exp), FUN = function(SAMPLE){
        # setTxtProgressBar(pb, SAMPLE)
        z <- obtain_max_gene_gMCS_sample(gene.exp.gMCSs[,SAMPLE], gMCSs.ENSEMBL.txt.list.num, genes.gMCSs.ENSEMBL, obtain_max_gene_gMCS)
        return(z)
      })
      # close(pb)
      gene.first.ENSEMBL <- do.call(cbind, gene.first.ENSEMBL)
      
      on.exit(plan(oplan), add = TRUE)
    }
    
    # once we now the gene, we can pass it to the sample
    for (SAMPLE in 1:ncol(gene.exp)){
      gene.first.exp[,SAMPLE] <- gene.exp.gMCSs[gene.first.ENSEMBL[,SAMPLE],SAMPLE]
    }
    colnames(gene.first.exp) = colnames(gene.first.ENSEMBL) =
      colnames(gene.exp)
    
    # Calculate threshold based on quantile ###
    ratio_threshold <- rep(NA,dim(gene.first.exp)[2])
    names(ratio_threshold) <- colnames(gene.first.exp)
    
    for (i in 1:dim(gene.first.exp)[2]){
      
      aux <- which(!duplicated(gene.first.ENSEMBL[,i]))
      aux <- gene.first.exp[aux,i]
      
      ratio_threshold[i] <- quantile(aux[aux>threshold_exp], gmcsTH_perc)
    }
    
    # define the gene ratio over the threshold ###
    gene.ratio <- gene.exp
    
    if (dim(gene.ratio)[1]>0) {
      for (i in 1:dim(gene.ratio)[2]) {
        gene.ratio[,i] <- gene.exp[,i]/ratio_threshold[i]
      }
    }
    gene.ON.OFF = log2(gene.ratio) >= threshold_logFC
    
    # Export the results
    return(list(
      localT2 = NULL,
      ratio_threshold = ratio_threshold,
      gene.ratio = gene.ratio,
      gene.ON.OFF = gene.ON.OFF
    ))
    
  } else if (thresholding_methodology == "singleTH"){
    ###########################  singleTH methodology   ###########################
    
    
    
    print("Perform the singleTH thesholding")
    
    # define the gene ratio over the threshold ###
    gene.ratio <- gene.exp/singleTH
    
    gene.ON.OFF = log2(gene.ratio) >= threshold_logFC
    
    # Export the results
    return(list(
      localT2 = NULL,
      gene.ratio = gene.ratio,
      gene.ON.OFF = gene.ON.OFF
    ))
    
  } else if (thresholding_methodology == "localT2"){
    ###########################  localT2 methodology   ###########################
    
    # Calculate localT2 thresholds  ###
    print("Perform the localT2 thesholding")
    
    localT2 <- list()
    
    for (cm in levels(sample.cohort)){
      print(paste(cm))
      
      gene.cohort.TPM <- gene.exp[genes.in.set, sample.cohort==cm, drop = F]
      
      localT2[[cm]] <- list()
      localT2[[cm]]$OFF <- quantile(gene.cohort.TPM, 0.25, na.rm = T)
      localT2[[cm]]$ON <- quantile(gene.cohort.TPM, 0.75, na.rm = T)
      localT2[[cm]]$MAYBEON <- apply(gene.cohort.TPM,1,mean)
      
      # print(localT2[[cm]][c(1,2)])
      # print(localT2[[cm]]$MAYBEON[1:5])
    }
    
    
    # Define genes in cohort as ON or OFF ###
    gene.ON.OFF <- gene.exp[genes.in.set, , drop = F]*0
    
    for (cm in levels(sample.cohort)){
      print(paste(cm))
      gene.aux.TPM <- gene.exp[genes.in.set, sample.cohort==cm, drop = F]
      gene.aux.ON.OFF <-  gene.aux.TPM*0
      
      gene.aux.TH.MAYBEON <- gene.aux.TPM*0
      for (i in 1:dim(gene.aux.TH.MAYBEON)[2])
      { gene.aux.TH.MAYBEON[,i] <- localT2[[cm]]$MAYBEON[rownames(gene.aux.ON.OFF)] }
      
      gene.aux.ON.OFF <- (gene.aux.TPM >= gene.aux.TH.MAYBEON)
      gene.aux.ON.OFF[gene.aux.TPM > localT2[[cm]]$ON] <- TRUE
      gene.aux.ON.OFF[gene.aux.TPM < localT2[[cm]]$OFF] <- FALSE
      
      gene.ON.OFF[, sample.cohort==cm] <-  gene.aux.ON.OFF
      rm(list = c("gene.aux.TPM", "gene.aux.ON.OFF"))
    }
    
    # Export the results
    return(list(
      localT2 = localT2,
      gene.ratio = NULL,
      gene.ON.OFF = gene.ON.OFF
    ))
    
  } else if (thresholding_methodology == "fastcormics"){
    ###########################  fastcormics methodology   ###########################
    
    # Calculate lfastcormics thresholds  ####
    print("Perform the fastcormics thesholding")
    
    gene.ratio <- gene.exp*0
    
    # define the set of genes
    gene.exp.set <- gene.exp[genes.in.set,]
    rownames(gene.exp.set) <- genes.in.set
    
    pb <- txtProgressBar(max = ncol(gene.exp), style = 3)
    
    for (SAMPLE in 1:ncol(gene.exp)){
      setTxtProgressBar(pb, SAMPLE)
      
      # browser();
      # print(SAMPLE)
      
      data.mat <- gene.exp.set[,SAMPLE]
      data.mat <- as.matrix(data.mat)
      data.mat <- data.mat[!is.na(data.mat)]
      if (log2.fastcormics){
        data.mat[data.mat<1e-6] <- 1e-6; data.mat <- log2(data.mat)
      }
      
      
      # integration step
      step = 1e-3
      limits = c(min(data.mat)-2,max(data.mat)+2)  # change limits to include all the values
      limits = round(limits, ceiling(-log10(step)))
      
      x <- seq(limits[1],limits[2],step)
      x <- round(x, ceiling(-log10(step)))
      
      # fit to density plot
      fit <- density(data.mat)
      
      # the dataframe will store everything
      df <- data.frame(x = x, hibrid = approx(fit$x, fit$y, x)$y)
      
      # calculate right hybrid
      fit2.mu <- df$x[which.max(df$hibrid)]
      fit2 <- df[df$x > fit2.mu,] # select half gaussian
      fit2$x <- fit2$x - fit2.mu
      fit2 <- data.frame(x = c(-rev(fit2$x[-1]),fit2$x) + fit2.mu, # create entire gaussian
                         hibrid_right = c(rev(fit2$hibrid[-1]), fit2$hibrid))
      
      df$x <- round(df$x, ceiling(-log10(step)))
      fit2$x <- round(fit2$x, ceiling(-log10(step)))
      df <- merge(df, fit2)
      df[is.na(df)] <- 0
      
      df$hibrid_left <- df$hibrid - df$hibrid_right
      
      df$hibrid[df$hibrid < 1e-3] <- 0
      df$hibrid_right[df$hibrid_right < 1e-3] <- 0
      df$hibrid_left[df$hibrid_left < 1e-3] <- 0
      
      # estimate the values
      xx <- df$x
      yy <- df$hibrid_right/sum(df$hibrid_right)*length(df$hibrid_right)
      zz <-  c(0.5, sum(xx*yy)/sum(yy), sqrt(sum((xx-sum(xx*yy)/sum(yy))^2 * yy)/sum(yy))) # initial estimations
      result = optim(zz,    function(pars){
        pred_pdf = pars[1]*dnorm(xx, mean = pars[2], sd = pars[3])
        err = sum((df$hibrid_right - pred_pdf)^2)
      }, method = "Nelder-Mead")
      df$calculated_right <- result$par[1] * dnorm(x = df$x, mean = result$par[2], sd = result$par[3])
      
      
      xx <- df$x
      yy <- df$hibrid_left/sum(df$hibrid_left)*length(df$hibrid_left)
      zz <-  c(0.5, sum(xx*yy)/sum(yy), sqrt(sum((xx-sum(xx*yy)/sum(yy))^2 * yy)/sum(yy))) # initial estimations
      result = optim(zz,    function(pars){
        pred_pdf = pars[1]*dnorm(xx, mean = pars[2], sd = pars[3])
        err = sum((df$hibrid_left - pred_pdf)^2)
      }, method = "Nelder-Mead")
      df$calculated_left <- result$par[1] * dnorm(x = df$x, mean = result$par[2], sd = result$par[3])
      
      
      df$hibrid[df$hibrid < 1e-3] <- 0
      df$hibrid_right[df$hibrid_right < 1e-3] <- 0
      df$hibrid_left[df$hibrid_left < 1e-3] <- 0
      df$calculated_right[df$calculated_right < 1e-3] <- 0
      df$calculated_left[df$calculated_left < 1e-3] <- 0
      
      df <- df[apply(df[,-1],1,sum)>1e-3,]
      
      if (plot.fastcormics){
        p <- ggplot(df) +
          geom_histogram(aes(x = x, y = ..density..),
                         data = data.frame(x = data.mat), bins = 100, alpha = 0.5, col = "grey20",fill = "grey80") +
          geom_line(aes(x = x, y = hibrid), size = 1.5, linetype = "dashed", color = "black") +
          geom_line(aes(x = x, y = hibrid_right), size = 1.25, linetype = "dashed", color = "red") +
          geom_line(aes(x = x, y = hibrid_left), size = 1.25, linetype = "dashed", color = "blue") +
          geom_line(aes(x = x, y = calculated_right), size = 1, linetype = "solid", color = "red") +
          geom_line(aes(x = x, y = calculated_left), size = 1, linetype = "solid", color = "blue") +
          theme_bw() + theme(legend.position="top") + xlab("log2(exp)") +
          ggtitle(paste0("Density plot and fitted curves of "), colnames(gene.exp)[SAMPLE])
        ggsave(paste0(file.path(plot.fastcormics.dir,colnames(gene.exp)[SAMPLE]),'.png'), scale = 1.5, plot = p)
      }
      
      mu1 <- df$x[which.max(df$calculated_left)]
      mu2 <- df$x[which.max(df$calculated_right)]
      aux <- gene.exp[,SAMPLE]; if (log2.fastcormics) {aux <- log2(aux)}
      aux = ((aux - mu1)/(mu2-mu1))*2-1
      aux[is.na(aux)] <- 0
      # aux[aux >= (+1)] <- (+1)
      # aux[aux <= (-1)] <- (-1)
      gene.ratio[,SAMPLE] <- aux
      
    }
    close(pb)
    
    
    gene.ON.OFF = gene.ratio >= 1
    
    # Export the results
    return(list(
      localT2 = NULL,
      gene.ratio = gene.ratio, # -1 if less than the mean of the left, +1
      gene.ON.OFF = gene.ON.OFF
    ))
    
  } else {stop("no thresholding methodology defined")}
  
}



