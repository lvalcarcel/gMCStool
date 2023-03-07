# rm(list = ls())
# gc()
# 
# setwd("C:/Users/lvalcarcel/OneDrive - Tecnun/Postdoc/002 - MM Human-GEM NY - CoMMpass - CCLE EBI 2021-06/00 - GENERATE APP/gMCStool-shinyApp-v5/CODE_TEST_HEATMAP")
# 
# library(ggpubr)
# 
# values = list(readRDS("gMCStool_heatmap_data_2023_03_06_13h09m26s.RDS"),
#               readRDS("gMCStool_heatmap_data_2023_03_06_13h09m49s.RDS"))
# 
# values = values[[1]]
# 
# gene.exp = values$gene.exp
# gene.ratio = values$gene.ratio
# localT2 = values$localT2
# gMCS.info.all = values$gMCS.info.all
# gMCS.info.target = values$gMCS.info.target
# sample.class = values$sample.class
# sample.class.target = values$sample.class.target
# sample.cohort = values$sample.cohort
# flag_include_boxplots_heatmap = values$flag_include_boxplots_heatmap
# flag_show_SYMBOL = values$flag_show_SYMBOL
# colors_heatmap = values$colors_heatmap
# colors_annotation = values$colors_annotation
# hline_value = 1
# TABLE_GMCS_MODE = values$TABLE_GMCS_MODE


ShowHeatmap <- function(gene.exp = NULL,
                        gene.ratio = NULL,
                        localT2 = NULL,
                        gMCS.info.all = NULL,
                        gMCS.info.target = NULL,
                        sample.class = NULL,
                        sample.class.target = NULL,
                        sample.cohort = NULL,
                        flag_include_boxplots_heatmap = T,
                        flag_show_SYMBOL = T,
                        colors_heatmap = "#902EBB_white_#F4831B",
                        colors_annotation = 1,
                        hline_value = 1,
                        # isDoubleKO = F,
                        TABLE_GMCS_MODE){
  
  # browser()
  
  colors_annotation <- list(rainbow, scales::hue_pal(), scales::hue_pal())[[as.numeric(as.character(colors_annotation))]]
  sample.class.target <- sample.class.target[1:length(levels(sample.class))]
  sample.class.target[is.na(sample.class.target)] <- F
  
  # transform gene expression into gene.ratio if localT2 is provided
  if (!is.null(localT2)){
    # extrat the genes and produce the log2ratio
    gene.set <- gene.exp[names(localT2[[1]]$MAYBEON),]*0
    for (k in 1:ncol(gene.exp))
    {
      aux <- localT2[[sample.cohort[k]]]$MAYBEON
      
      aux2 <- log2(gene.exp[names(localT2[[1]]$MAYBEON),k]/aux)
      aux2[gene.exp[names(localT2[[1]]$MAYBEON),k] > localT2[[sample.cohort[k]]]$ON] <- (+Inf)
      aux2[gene.exp[names(localT2[[1]]$MAYBEON),k] < localT2[[sample.cohort[k]]]$OFF] <- (-Inf)
      gene.set[,k] <- aux2
    }
  } else {
    gene.set = log2(gene.ratio)
  }
  
  # extract the neccesary info from the target
  if ("ENSEMBL.1" %in% colnames(gMCS.info.target)){
    gen.ENSEMBL <- c(gMCS.info.target$ENSEMBL.1, gMCS.info.target$ENSEMBL.2)
    gen.SYMBOL <- c(gMCS.info.target$SYMBOL.1, gMCS.info.target$SYMBOL.2)
    index <- as.numeric(as.character(gMCS.info.target$gMCS))
  } else {
    gen.ENSEMBL <- gMCS.info.target$ENSEMBL
    gen.SYMBOL <- gMCS.info.target$SYMBOL
    index <- as.numeric(as.character(gMCS.info.target$gMCS))
  }
  COLNAMES <- colnames(gene.exp)
  
  # gmcs <- unique(as.character(gMCSs.ENSEMBL[[n]][index,]))
  gmcs.ENSEMBL <- unique(unlist(strsplit(as.character(gMCS.info.all$gMCSs.ENSEMBL.txt[index]),'--')))
  gmcs.ENSEMBL <- gmcs.ENSEMBL[!is.na(gmcs.ENSEMBL)]
  gene.set <- (as.matrix(gene.set[gmcs.ENSEMBL,]))
  gene.set.exp <- as.matrix(gene.exp[gmcs.ENSEMBL,])
  # colnames(gene.set) <- sample.class
  rownames(gene.set) <- gmcs.ENSEMBL
  # try({
  if (length(gmcs.ENSEMBL)==2 & length(gen.ENSEMBL)==1) { 
    gene.set <- gene.set[c(gen.ENSEMBL, gmcs.ENSEMBL[!(gmcs.ENSEMBL %in% gen.ENSEMBL)]),]
  } else if (length(gmcs.ENSEMBL)==3 & length(gen.ENSEMBL)==2) { 
    gene.set <- gene.set[c(gen.ENSEMBL, gmcs.ENSEMBL[!(gmcs.ENSEMBL %in% gen.ENSEMBL)]),]
  } else if (length(gmcs.ENSEMBL)>3) {
    gene.set2 <- gene.set[!(gmcs.ENSEMBL %in% gen.ENSEMBL),]
    gene.set2[gene.set2>10] <- 10
    gene.set2[gene.set2<(-10)] <- (-10)
    gene.set <- gene.set[c(gen.ENSEMBL,names(sort(apply(gene.set2,1,mean, na.rm = T),decreasing = T))),]
  }
  # })
  gmcs.SYMBOL  <- gMCS.info.all$table.genes.HumanGEM$SYMBOL[match(rownames(gene.set),gMCS.info.all$table.genes.HumanGEM$ENSEMBL )]
  gmcs.ENSEMBL <- rownames(gene.set)
  
  
  if (length(gmcs.ENSEMBL)==1)  {
    gene.set <- matrix(gene.ratio[gmcs.ENSEMBL,], nrow = 1)
    gene.set.exp <- matrix(gene.exp[gmcs.ENSEMBL,], nrow = 1)
    colnames(gene.set) <- COLNAMES
    colnames(gene.set.exp) <- COLNAMES
  } else {
    gene.set <- as.matrix(gene.set)
    gene.set.exp <- as.matrix(gene.set.exp[rownames(gene.set),])
  }
  
  
  # Separate classes if there some of the classes are target
  sample.class.2 <- as.character(sample.class)
  
  for (cl in levels(sample.class)[sample.class.target>0.5]) {
    cl.responder <- paste0(cl, ".responder")
    cl.non.responder <- paste0(cl, ".non.responder")
    sample.class.2[sample.class==cl] <- cl.non.responder
    if (length(gen.ENSEMBL) == 1){
      gene.set.rest <- gene.set[-1,sample.class==cl, drop = F]
      # if (length(gmcs.ENSEMBL)==2)  {  gene.set.rest <- matrix(gene.set.rest, nrow = 1) }
      sample.class.2[sample.class==cl][gene.set[1,sample.class==cl]>0 & apply(gene.set.rest,2,max)<0] <- cl.responder
    } else {
      gene.set.rest <- gene.set[-c(1,2),sample.class==cl, drop = F]
      sample.class.2[sample.class==cl][gene.set[1,sample.class==cl]>0 & gene.set[2,sample.class==cl]>0 & apply(gene.set.rest,2,max)<0] <- cl.responder
    }
    
  }
  
  sample.class.2.levels <- c()
  for(i in 1:length(levels(sample.class))){
    cl <- levels(sample.class)[i]
    if (sample.class.target[i]>0.5){
      sample.class.2.levels <- c(sample.class.2.levels, paste0(cl, ".non.responder"), paste0(cl, ".responder"))
    } else {
      sample.class.2.levels <- c(sample.class.2.levels, cl)
    }
  }
  # sample.class.2.levels <- sample.class.2.levels[sample.class.2.levels %in% sample.class.2]
  
  sample.class.2 <- factor(sample.class.2, levels = sample.class.2.levels)
  
  idx <- order(sample.class.2)
  sample.class <- sample.class[idx]
  sample.class.2 <- sample.class.2[idx]
  sample.cohort <- sample.cohort[idx]
  gene.set <- gene.set[,idx]
  gene.set.exp <- gene.set.exp[,idx]
  COLNAMES <- COLNAMES[idx]
  
  if (length(gmcs.ENSEMBL)==1)  {
    gene.set <- matrix(gene.set, nrow = 1)
    gene.set.exp <- matrix(gene.set.exp, nrow = 1)
    colnames(gene.set) <- COLNAMES
    colnames(gene.set.exp) <- COLNAMES
  }
  
  
  # set rownames 
  if (flag_show_SYMBOL) { 
    rownames(gene.set) <- gmcs.SYMBOL
    rownames(gene.set.exp) <- gmcs.SYMBOL
    # gmcs <- gmcs.SYMBOL
  } else {
    rownames(gene.set) <- gmcs.ENSEMBL
    rownames(gene.set.exp) <- gmcs.ENSEMBL
    # gmcs <- gmcs.ENSEMBL
  }
  
  gene.set.list <- lapply(levels(sample.cohort), function(x){gene.set[,sample.cohort == x]})
  names(gene.set.list) <- levels(sample.cohort)
  
  colnames_aux <- lapply(levels(sample.cohort), function(x){colnames(gene.set)[sample.cohort == x]})
  names(colnames_aux) <- levels(sample.cohort)
  
  gene.set.exp.list <- lapply(levels(sample.cohort), function(x){gene.set.exp[,sample.cohort == x]})
  names(gene.set.exp.list) <- levels(sample.cohort)
  
  
  
  if (length(gmcs.ENSEMBL)==1){
    for (ch in levels(sample.cohort)) {
      # colnames_aux <- colnames(gene.set.list[[ch]])
      gene.set.list[[ch]] <- matrix(unlist(gene.set.list[[ch]]), nrow = 1)
      gene.set.exp.list[[ch]] <- matrix(unlist(gene.set.exp.list[[ch]]), nrow = 1)
      colnames(gene.set.list[[ch]]) <- colnames_aux[[ch]]
      colnames(gene.set.exp.list[[ch]]) <- colnames_aux[[ch]]
      if (flag_show_SYMBOL) { rownames(gene.set.list[[ch]]) <- gmcs.SYMBOL } else {rownames(gene.set.list[[ch]]) <- gmcs.ENSEMBL}
      if (flag_show_SYMBOL) { rownames(gene.set.exp.list[[ch]]) <- gmcs.SYMBOL } else {rownames(gene.set.exp.list[[ch]]) <- gmcs.ENSEMBL}
    }
  }
  
  
  
  sample.class.list <- lapply(levels(sample.cohort), function(ch){
    z <- sample.class[sample.cohort == ch]
    z <- factor(as.character(z),levels=unique(as.character(z), fromLast = T))
    return(z)
  })
  names(sample.class.list) <- levels(sample.cohort)
  
  
  sample.class.2.list <- lapply(levels(sample.cohort), function(ch){
    z <- as.character(sample.class.2[sample.cohort == ch])
    z <- factor(z,levels=levels(sample.class.2)[levels(sample.class.2) %in% z])
    return(z)
  })
  names(sample.class.2.list) <- levels(sample.cohort)
  
  
  
  # Generate colnames annotation (prepare to have same length)
  
  ANNOTATION_TEXT <- list()
  for (ch in levels(sample.cohort)) {
    # change column names to include the sample class
    if (ncol(gene.set.list[[ch]])>20){
      ANNOTATION_TEXT[[ch]] <- rep("",length(sample.class.2.list[[ch]]))
      ANNOTATION_TEXT[[ch]][round(unlist(lapply(as.list(levels(sample.class.2.list[[ch]])),function(x){median(which(sample.class.2.list[[ch]][order(sample.class.2.list[[ch]])]==x))})))] <- levels(sample.class.2.list[[ch]])
    } else {
      ANNOTATION_TEXT[[ch]] <- colnames(gene.set.list[[ch]])
    }
  }
  
  nn <- max(unlist(lapply(ANNOTATION_TEXT, nchar)))
  ANNOTATION_TEXT <- lapply(ANNOTATION_TEXT, function(aux) {
    # aux <- ANNOTATION_TEXT$CCLE
    for(i in 1:length(aux)) { aux[i] <- paste0(aux[i], paste0(rep(" ", ceiling((nn-nchar(aux[i]))*1.25)), collapse = ""))}
    return(aux)
  })
  
  
  ## Generate a Heatmap for each cohort
  plot.list <- list()
  
  
  for (ch in levels(sample.cohort)) {
    # ch <- levels(sample.cohort)[1]
    # define the colors
    COLORS <- colors_annotation(length(levels(sample.class.2)))
    names(COLORS) <- levels(sample.class.2)
    
    if (length(gmcs.ENSEMBL)>1){
      # order by gene essentiality
      
      if (length(gen.ENSEMBL) == 1){
        gene.set.rest <- gene.set.list[[ch]][-1,,drop = F]
        idx <- order(gene.set.list[[ch]][1,] - apply(gene.set.rest,2,max), decreasing = F)
      } else {
        gene.set.rest <- gene.set.list[[ch]][-c(1:2),,drop = F]
        idx <- order(apply(gene.set.list[[ch]][1:2,,drop = F],2,min) - apply(gene.set.rest,2,max), decreasing = F)
      }
      # gene.set.rest <- gene.set.list[[ch]][-1,,drop = F]
      # idx <- order(gene.set.list[[ch]][1,] - apply(gene.set.rest,2,max), decreasing = F)
      sample.class.list[[ch]] <- sample.class.list[[ch]][idx]
      sample.class.2.list[[ch]] <- sample.class.2.list[[ch]][idx]
      gene.set.list[[ch]] <- gene.set.list[[ch]][,idx]
      
      # order by class
      idx <- order(sample.class.2.list[[ch]])
      sample.class.list[[ch]] <- sample.class.list[[ch]][idx]
      sample.class.2.list[[ch]] <- sample.class.2.list[[ch]][idx]
      gene.set.list[[ch]] <- gene.set.list[[ch]][,idx]
    } else {
      # order by class
      idx <- order(sample.class.2.list[[ch]])
      COLNAMES_aux <- colnames(gene.set.list[[ch]])
      sample.class.list[[ch]] <- sample.class.list[[ch]][idx]
      sample.class.2.list[[ch]] <- sample.class.2.list[[ch]][idx]
      gene.set.list[[ch]] <- matrix(gene.set.list[[ch]][,idx], nrow = 1)
      colnames(gene.set.list[[ch]]) <- COLNAMES_aux
      # set rownames 
      if (flag_show_SYMBOL) { rownames(gene.set.list[[ch]]) <- gmcs.SYMBOL } else {rownames(gene.set.list[[ch]]) <- gmcs.ENSEMBL}
    }
    
    ha <- data.frame(class = sample.class.2.list[[ch]], row.names = colnames(gene.set.list[[ch]]))
    
    
    if (sum(sample.class.target) > 0.5) { 
      aux <- intersect(levels(sample.class)[sample.class.target>0.5], sample.class.list[[ch]])
      if (TABLE_GMCS_MODE=="percentage") {
        aux2 <- gMCS.info.target[,aux]
      } else {
        aux2 <- (gMCS.info.target[,levels(sample.class)]/table(sample.class))[sample.class.target>0.5]
      }
      if (length(aux)>1){
        aux <- paste(aux, " ", round(aux2*100, 1), "%", sep = "", collapse = ", ")
      } else {
        aux <- paste0(round(aux2*100, 1), "%")
      }
      aux <- paste0(" (", aux, ")")
    } else {aux <- ""}
    
    plot.list[[ch]] <- pheatmap(gene.set.list[[ch]],
                                cluster_rows = F, cluster_cols = F,
                                gaps_col = cumsum(table(sample.class.list[[ch]])),
                                show_rownames = TRUE, show_colnames = T,  
                                labels_col = ANNOTATION_TEXT[[ch]], 
                                color = colorRampPalette(unlist(strsplit(colors_heatmap,'_')))(100), 
                                breaks = seq(-1,+1, length.out = 100),
                                annotation_col = ha,
                                legend = ch==levels(sample.cohort)[length(levels(sample.cohort))],
                                annotation_legend = ch==levels(sample.cohort)[length(levels(sample.cohort))],
                                annotation_colors = list(class = COLORS),
                                border_color = ifelse(ncol(as.matrix(gene.set.list[[ch]]))<=20,"grey60",NA),
                                main = paste0(ch, aux), 
                                scale = "none",
                                silent = T)[[4]]
  }
  
  pp_heatmaps_widths <- log10(table(sample.cohort))+2
  pp_heatmaps_widths[length(pp_heatmaps_widths)] <- pp_heatmaps_widths[length(pp_heatmaps_widths)] +1
  
  pp_heatmaps <- ggpubr::ggarrange(plotlist = plot.list, nrow = 1, ncol = length(levels(sample.cohort)),
                                   widths = pp_heatmaps_widths)
  
  
  
  ## Generate a boxplot for each cohort
  
  if (flag_include_boxplots_heatmap){
    plot.list.2 <- list()
    levels.list.2 <- list()
    
    for (ch in levels(sample.cohort)) {
      # ch <- levels(sample.cohort)[1]
      # define the colors
      COLORS <- colors_annotation(length(levels(sample.class.2.list[[ch]])))
      names(COLORS) <- levels(sample.class.2.list[[ch]])
      
      COLORS <- colors_annotation(length(levels(sample.class.2)))
      names(COLORS) <- levels(sample.class.2)
      
      if (length(gmcs.ENSEMBL)==1)  {
        sdf <- data.frame(Var1 = rownames(gene.set.exp.list[[ch]]),
                          Var2 = colnames(gene.set.exp.list[[ch]]),
                          value = as.numeric(as.matrix(gene.set.exp.list[[ch]])))
      } else {
        sdf <- reshape2::melt(as.matrix(gene.set.exp.list[[ch]]))
      }
      sdf$Var3 <- match(sdf$Var2,colnames(gene.set.exp))
      sdf$class <- sample.class.2[sdf$Var3]
      sdf$exp <- sdf$value
      sdf$log2exp <- log2(sdf$exp+1)
      
      levels.list.2[[ch]] <- length(unique(sdf$class))
      
      plot.list.2[[ch]] <- ggboxplot(sdf, x = "Var1", y  ="exp", fill = "class") + theme_bw() + 
        xlab("") + ylab("expression") + theme(legend.position = "bottom") + 
        # scale_y_continuous(trans = "log1p", limits = c(0,NA))+
        scale_fill_manual(values = COLORS, breaks = names(COLORS)) + ggtitle(ch) + 
        theme(legend.title = element_blank())
      
      if (!is.null(localT2)){
        
        dLines <- data.frame(X = 1:length(gmcs.ENSEMBL)-0.45,
                             Y = localT2[[ch]]$MAYBEON[gmcs.ENSEMBL],
                             Xend = 1:length(gmcs.ENSEMBL)+0.45,
                             Yend=localT2[[ch]]$MAYBEON[gmcs.ENSEMBL],
                             group = rep("MAYBE ON", length(gmcs.ENSEMBL)))
        
        plot.list.2[[ch]] <- plot.list.2[[ch]] + 
          geom_hline(yintercept = c(localT2[[ch]]$ON, localT2[[ch]]$OFF),  col = "red", linetype = "dotted", size = 1) + 
          geom_segment(data = dLines, aes(x = X, xend = Xend, y = Y, yend = Yend),color="red",size=1,linetype="twodash", show.legend = F)
      } else if (!is.null(hline_value) & !is.na(hline_value)){
        plot.list.2[[ch]] <- plot.list.2[[ch]] + 
          geom_hline(yintercept = hline_value, col = "red") + guides(fill = guide_legend(nrow = 1)) 
      }
      
      # detect if data us log2 transformed
      qx <- as.numeric(quantile(gene.exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
      LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
      
      # Add p-values of `stat.test` and `stat.test2`
      stat.test <- sdf
      if (LogC) {stat.test$exp <- log2(sdf$exp+1)} #else {stat.test$log2exp <- sdf$exp}
      if (min(table(sdf$class)[table(sdf$class)>0])/length(gmcs.ENSEMBL)>1){
        try({
          stat.test <- stat.test %>%
            group_by(Var1) %>%
            t_test(log2exp~class) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj")
          # 1. Add stat.test
          stat.test <- stat.test %>%
            add_xy_position(x = "Var1", dodge = 0.8)
          # if (ch == "Bcell and MM") {
          # 2. reduce the p-values to only responders
          # browser()
          stat.test.2 <- stat.test
          stat.test <- stat.test %>% filter(grepl("responder", group1))  %>% filter(grepl("responder", group2))
          stat.test <- stat.test %>% filter(grepl("responder", group1))  %>% filter(grepl("responder", group2))
          stat.test <- stat.test %>% filter(sub(".responder", "", sub(".non.responder", "", group1)) == sub(".responder", "", sub(".non.responder", "", group2)))
          stat.test$y.position <- stat.test.2$y.position[unlist(lapply(unique(stat.test$Var1), function(x) which(stat.test.2$Var1==x)[1:table(stat.test$Var1)[x]]))]
          # max(stat.test %>% filter(!(grepl("responder", group1) & grepl("responder", group2))) %>% select(y.position))
          # }
          # Add p-values onto the bar plots
          # stat.test <- stat.test %>%
          #   add_xy_position(fun = "mean_sd", x = "Var1", dodge = 0.8) 
          plot.list.2[[ch]] <- plot.list.2[[ch]] + stat_pvalue_manual( stat.test,  label = "p.adj.signif",  tip.length = 0 )
        })
      }
      
      if (!is.null(localT2)){
        second_axis_aux <- sec_axis(~.*1, name="",
                                    breaks = c(localT2[[ch]]$ON, localT2[[ch]]$OFF),
                                    labels = c("ON", "OFF"))
        if (LogC) { 
          plot.list.2[[ch]] <- plot.list.2[[ch]] + scale_y_continuous(trans = "log1p", limits = c(0,NA), sec.axis = second_axis_aux)
        } else {
          plot.list.2[[ch]] <- plot.list.2[[ch]] + scale_y_continuous(sec.axis = second_axis_aux) 
        }
      } else {
        if (LogC) { 
          plot.list.2[[ch]] <- plot.list.2[[ch]] + scale_y_continuous(trans = "log1p", limits = c(0,NA))
        }
      }
      
    }
    pp_boxplot <-  ggpubr::ggarrange(plotlist = plot.list.2, nrow = 1, ncol = length(levels(sample.cohort)),
                                     widths = log10(unlist(levels.list.2))+2, common.legend = T, legend = "bottom")
  }
  
  
  
  if (flag_include_boxplots_heatmap){
    return(list(heatmap = pp_heatmaps, boxplot = pp_boxplot))
    # ggarrange(pp_heatmaps, pp_boxplot, ncol = 1)
  } else {
    return(list(heatmap = pp_heatmaps))
    # ggarrange(pp_heatmaps)
  }
}
