ShowDotplot_DepMap <- function(DepMap.info.all, gMCS.info.all, gmcs_database, gene.target.info, by_gMCS = F,
                               database, database_filter_mode, database_filter_selected, threshold_value,
                               flag_database_filter_show_only, flag_LinearRegression = T, flag_show_SYMBOL = T,
                               flag_color_by_gMCS = F, database_unit){
  
  # browser()
  
  # remove empty database
  # database_filter_selected <- database_filter_selected[database_filter_selected != "---"]
  
  # extract the neccesary info from the target
  gene.ENSEMBL <- gene.target.info$ENSEMBL
  gene.SYMBOL <- gene.target.info$SYMBOL
  
  
  # test the data
  if(gene.ENSEMBL %in% gMCS.info.all$gMCSs.ENSEMBL.txt) { return(tibble(x = 1, y = 1) %>% ggplot(aes(x, y)) + theme_void() + geom_text(label ="This is an essential gene, it has no partner")) }
  # if(sum(gMCS.info.all$gMCSs.ENSEMBL.max[,gene.ENSEMBL])==1){by_gMCS = F; flag_color_by_gMCS = F}
  
  if (by_gMCS) {
    gmcs.ENSEMBL <- strsplit(gMCS.info.all$gMCSs.ENSEMBL.txt[as.numeric(as.character(gene.target.info$gMCS))],'--')[[1]]
    gmcs.ENSEMBL <- setdiff(gmcs.ENSEMBL, gene.ENSEMBL)
    gmcs.SYMBOL <- strsplit(gMCS.info.all$gMCSs.SYMBOL.txt[as.numeric(as.character(gene.target.info$gMCS))],'--')[[1]]
    gmcs.SYMBOL <- setdiff(gmcs.SYMBOL, gene.SYMBOL)
    
    sdf <- reshape2::dcast(DepMap.info.all$DepMapGeneExpression %>% filter(UNIT==database_unit), ENSEMBL~DepMap_ID, value.var = "logTPM", fun.aggregate = sum)
    sdf <- sdf %>% filter(ENSEMBL %in% gmcs.ENSEMBL)
    sdf <- reshape2::melt(apply(sdf[,-1],2,max))
    sdf <- data.frame(DepMap_ID = rownames(sdf),
                      logTPM = sdf$value,
                      ENSEMBL = gene.ENSEMBL,
                      UNIT = database_unit)
    
  } else {
    sdf <- DepMap.info.all$DepMapExpressionByGene %>% 
      filter(ENSEMBL==gene.ENSEMBL) %>% 
      filter(gMCS.database==gmcs_database) %>% 
      filter(UNIT==database_unit) 
  }
  
  
  if (flag_color_by_gMCS & !by_gMCS) {
    
    aux.fun <-  function(x){
      if (length(x)==1){ 
        paste0("",x , "")
      } else {
        paste0("max(", paste(x, collapse = ", "), ")")
      }}
    
    # define the gMCSs in which this gene participates
    idx <- sdf$gmcs.idx
    if (flag_show_SYMBOL) {
      # browser()
      gMCSs.SYMBOL <- gMCS.info.all$gMCSs.SYMBOL.txt[idx]
      gMCSs.SYMBOL <- lapply(gMCSs.SYMBOL, function(x){unlist(strsplit(x, "--"))})
      
      gMCSs.SYMBOL <- lapply(gMCSs.SYMBOL, function(x){x[x!=gene.SYMBOL]})
      sdf$gMCSs.to.show <- unlist(lapply(gMCSs.SYMBOL, aux.fun))
    } else {
      gMCSs.ENSEMBL <- gMCS.info.all$gMCSs.ENSEMBL.txt[idx]
      gMCSs.ENSEMBL <- lapply(gMCSs.ENSEMBL, function(x){unlist(strsplit(x, "--"))})
      
      gMCSs.ENSEMBL <- lapply(gMCSs.ENSEMBL, function(x){x[x!=gene.ENSEMBL]})
      
      sdf$gMCSs.to.show <- unlist(lapply(gMCSs.ENSEMBL, aux.fun))
    }
  } else {
    sdf$gMCSs.to.show <- "noColor"
  }
  
  # head(sdf)
  
  
  sdf.2 <- DepMap.info.all$DepMapEssentialityByGene %>% 
    filter(ENSEMBL==gene.ENSEMBL) %>% 
    filter(essentiality.database==database)
  
  sdf <- merge(sdf, sdf.2)
  sdf$logTPM <- as.numeric(as.character(sdf$logTPM))
  sdf$essentiality_score <- as.numeric(as.character(sdf$essentiality_score))
  
  # head(sdf)
  
  # add the filter selected
  sdf$isFilterSelected <- "rest"
  
  if (database_filter_mode!="none"){  # & database_filter_selected!="---"){
    sdf <- merge(sdf, DepMap.info.all$dictionary.CCLE[,c("DepMap_ID", database_filter_mode)])
    sdf$isFilterSelected[sdf[,database_filter_mode] %in% database_filter_selected] <- sdf[sdf[,database_filter_mode] %in% database_filter_selected, database_filter_mode]
  }
  sdf$isFilterSelected <- as.factor(sdf$isFilterSelected)
  
  
  XLAB <- paste0("Minimum expression of associated gMCSs partner genes \\[ ", sub("log2", "log_2", database_unit), " \\]")
  if (by_gMCS){
    if (flag_show_SYMBOL) { 
      XLAB <- paste0("Maximum expression of all partner genes: \\textit{", paste(gmcs.SYMBOL, collapse = ", "), "} \\[ ", sub("log2", "log_2", database_unit), " \\]")
    } else {
      XLAB <- paste0("Maximum expression of all partner genes: \\textit{", paste(gmcs.ENSEMBL, collapse = ", "), "} \\[ ", sub("log2", "log_2", database_unit), " \\]")
    }
  }
  
  
  # set the Y lab
  if (flag_show_SYMBOL) { 
    YLAB <- paste0("Essentiality of \\textit{",gene.SYMBOL , "} \\[ ",database," score \\]")
  } else {
    YLAB <- paste0("Essentiality of \\textit{",gene.ENSEMBL , "} \\[ ",database," score \\]")
  }
  
  ## Generate a Dotplot
  pp_dotplot <- ggplot(sdf, aes(x = logTPM, y = essentiality_score)) + 
    # scale_y_continuous() + 
    theme_classic() + 
    theme(text = element_text(size=12),
          axis.text.x = element_text(color='black', size=12),
          legend.title = element_blank(),
          # axis.ticks.x = element_blank(),
          axis.text.y = element_text(color='black', size=12),
          legend.text = element_text(color='black', size=12,  margin = margin(r = 30, unit = "pt")),
          # legend.spacing.x = unit(3, 'cm'),
          # legend.margin = margin(t = 0, r = 2, b = 1, l = 2, unit = "cm"),
          plot.title = element_text(color='black', size=20),
          # axis.line.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          strip.text = element_text(color='black', size=14), 
          axis.title = element_text(color = "black", size = 14, face = "bold"),
          legend.position = "none") + 
    xlab(latex2exp::TeX(XLAB)) + ylab(latex2exp::TeX(YLAB)) 
  
  if (!is.na(threshold_value)) {
    pp_dotplot <- pp_dotplot + geom_hline(yintercept = threshold_value, linetype = "dotdash", col = "red3", size = 1)
  }
  
  if (database_filter_mode!="none"){
    if (flag_database_filter_show_only){
      sdf3 <- sdf[sdf$isFilterSelected!="rest",]
      pp_dotplot <- pp_dotplot + 
        geom_point(mapping = aes(col = gMCSs.to.show, shape = isFilterSelected), size = 2, data = sdf3) +
        # scale_color_manual(values = c("grey50", "black"), breaks = c("0","1")) +
        scale_shape_manual(values = c(21, rep(16, length(database_filter_selected))), breaks = c("rest",database_filter_selected))  
      # scale_size_manual(values = c(1.5, 3), breaks = c("0","1"))
    } else {
      sdf3 <- sdf[sdf$isFilterSelected!="rest",]
      pp_dotplot <- pp_dotplot + 
        geom_point(mapping = aes(col = gMCSs.to.show, shape = isFilterSelected), size = 1.5, data = sdf) +
        geom_point(mapping = aes(col = gMCSs.to.show, shape = isFilterSelected), size = 3, data = sdf3) +
        # scale_color_manual(values = c("grey50", "black"), breaks = c("0","1")) +
        scale_shape_manual(values = c(21, rep(8, length(database_filter_selected))), breaks = c("rest",database_filter_selected))  
      # scale_size_manual(values = c(1.5, 3), breaks = c("0","1"))
    }
  } else {
    pp_dotplot <- pp_dotplot + geom_point(aes(col = gMCSs.to.show), size = 2)
  }
  
  if (flag_LinearRegression){
    if (flag_database_filter_show_only){
      sdf3 <- sdf[sdf$isFilterSelected!="rest",]
    } else { sdf3 <- sdf}
    pp_dotplot <- pp_dotplot +  geom_smooth(method = "lm", color = "black", fullrange = T,
                                            mapping = aes(x = logTPM, y = essentiality_score), formula = 'y ~ x',
                                            data = sdf3)
    
    z <- cor.test(sdf3$essentiality_score, sdf3$logTPM)
    pp_dotplot <- pp_dotplot + annotate("text", label = paste0("\u03C1 = ",round(z$estimate,3),"           \np-value = ",format.pval(z$p.value),"           "),
                                        x = Inf, y = -Inf, hjust = 1, vjust = -1, size = 4)
  }
  
  if (!flag_color_by_gMCS){
    pp_dotplot <- pp_dotplot + scale_color_manual(values = c("black"), breaks = c("noColor"))
  }
  
  # pp_dotplot
  return(pp_dotplot)
}
