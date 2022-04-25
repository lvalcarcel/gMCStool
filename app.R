## ==================================================================================== ##
# gMCStool Shiny App for predicting of essential genes.
# Copyright (C) 2021 Luis V. Valcarcel
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# You may contact the authors of this code, Luis V. Valcarcel, at <lvalcarcel@unav>
## ==================================================================================== ##
# Use saveRDS readRDS

gc()
cat("\f")

# Clean all object from workspace
rm(list = ls())
# browser()


# install.packages("rstudioapi")
library(rstudioapi) # Include it in helpers.R
setwd(dirname(getActiveDocumentContext()$path)) # Change into file directory

########################################################################################################
#  GLOBAL
########################################################################################################

source("helpers.R") # Here, all the libraries and scripts needed are loaded

if (TRUE){ # set to false if this is goinf to be run in AWS
  nWorkers = parallel::detectCores()-2  # Define maximum number of workers for parallel processing
  nBytesRAM = 4*1024*1024*1024
  nBytesRAM_OS = 8*1024*1024*1024
  nBytesRAM = max(nBytesRAM, benchmarkme::get_ram()-nBytesRAM_OS)  # Define maximum RAM for shiny
  num_max_samples = 1e6
} else {
  nWorkers = 2  # Define maximum number of workers for parallel processing
  nBytesRAM = 1.5*1000*1024*1024 # Define maximum RAM for shiny
  num_max_samples = 100
}

num_max_classes = 30 # define the maximum number of classes

RealTimeTables_mp4 = T
RealTimeTables_mp5 = F

# Load data  ####
gMCS.info.raw <- list()
gMCS.info.raw[["EssentialTasks_CultureMedium"]] <- new.env()
load("Data/gMCSs_EssentialTasks_CultureBiomass_combined_HumanGEMv1.4.0_ENSEMBL.rdata", envir = gMCS.info.raw[["EssentialTasks_CultureMedium"]])
gMCS.info.raw[["EssentialTasks_FullMedium"]] <- new.env()
load("Data/gMCSs_EssentialTasks_FullBiomass_combined_HumanGEMv1.4.0_ENSEMBL.rdata", envir = gMCS.info.raw[["EssentialTasks_FullMedium"]])
gMCS.info.raw[["Only_CultureMedium"]] <- new.env()
load("Data/gMCSs_CultureBiomass_combined_HumanGEMv1.4.0_ENSEMBL.rdata", envir = gMCS.info.raw[["Only_CultureMedium"]])
gMCS.info.raw[["Only_FullMedium"]] <- new.env()
load("Data/gMCSs_FullBiomass_combined_HumanGEMv1.4.0_ENSEMBL.rdata", envir = gMCS.info.raw[["Only_FullMedium"]])
gMCS.info.raw <- lapply(gMCS.info.raw, as.list)
# remove unnecessary field (gMCS.ENSEMBL.list)
gMCS.info.raw <- lapply(gMCS.info.raw, function(x){x[c("gMCSs.ENSEMBL.txt", "table.gMCSs", "gMCSs.ENSEMBL.length",
                                                       "gMCSs.ENSEMBL.mat", "genes.gMCSs.ENSEMBL", "table.genes.HumanGEM",
                                                       "gMCSs.ENSEMBL.txt.SYMBOL", "gMCSs.ENSEMBL")]})

# generate field with reduced information
gMCS.info.raw[["Custom_CultureMedium"]] <- gMCS.info.raw[["EssentialTasks_CultureMedium"]][c("table.gMCSs", "genes.gMCSs.ENSEMBL", "table.genes.HumanGEM")]
gMCS.info.raw[["Custom_FullMedium"]] <- gMCS.info.raw[["EssentialTasks_FullMedium"]][c("table.gMCSs", "genes.gMCSs.ENSEMBL", "table.genes.HumanGEM")]

gMCS.info.raw[["EssentialTasks_CultureMedium"]]$fullname <- "Essential Tasks and Growth on Ham's medium"
gMCS.info.raw[["EssentialTasks_FullMedium"]]$fullname <- "Essential Tasks and Growth on unconstrained medium"
gMCS.info.raw[["Only_CultureMedium"]]$fullname <- "Only growth on Ham's medium"
gMCS.info.raw[["Only_FullMedium"]]$fullname <- "Only growth on unconstrained medium"
gMCS.info.raw[["Custom_CultureMedium"]]$fullname <- "Selected metabolic tasks and growth on Ham's medium"
gMCS.info.raw[["Custom_FullMedium"]]$fullname <- "Selected metabolic tasks and growth on unconstrained medium"

gMCS.info.raw <- lapply(gMCS.info.raw, as.list)

metTasks <- readRDS("Data/metTasks_HumanGEMv1.4.0.RDS")

metsBiomass <- readRDS("Data/Metabolites_of_interest_gMCS_biomass_HumanGEMv1.4.0.RDS")



DepMap.info.all <- new.env()
load("Data/DepMap_info_genes_gMCS_HumanGEMv1.4.0.RData", envir = DepMap.info.all)
load("Data/DepMap_correlation_genes_HumanGEMv1.4.0.RData", envir = DepMap.info.all)
DepMap.info.all <- as.list(DepMap.info.all)
DepMap.info.all <- DepMap.info.all[setdiff(names(DepMap.info.all), "DepMapEssentiality")]
DepMap.info.all$DepMapGeneExpression <- reshape2::melt(DepMap.info.all$DepMapGeneExpression)
colnames(DepMap.info.all$DepMapGeneExpression) <- c("ENSEMBL", "SYMBOL","DepMap_ID", "logTPM", "UNIT")
DepMap.info.all$DepMapGeneExpression <- DepMap.info.all$DepMapGeneExpression %>% 
  mutate(ENSEMBL = as.factor(ENSEMBL)) %>% mutate(SYMBOL = as.factor(SYMBOL)) %>% 
  mutate(DepMap_ID = as.factor(DepMap_ID)) %>% mutate(UNIT = as.factor(UNIT))
DepMap.info.all <- lapply(DepMap.info.all, as.data.frame)



print(paste0("running gMCStool with ",nWorkers," cores and ",round(nBytesRAM/1024/1024/1024,1)," GBs of RAM"))

########################################################################################################
#  USER INTERFACE
########################################################################################################

# Options for Spinner
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)

# only activate to debug
# options(shiny.error = browser)

STYLE_NEXT_BUTTONS = "color: #ffffff; background-color: #77933C; border-color: #216b2c"
STYLE_NEXT_BUTTONS = "color: #ffffff; background-color: #e46c0a; border-color: #87430b"

# APPLICATION WEBPAGE ----------------------------
# Define UI for application that draws a histogram
ui <- tagList( 
  
  tags$head(includeScript("./www/google-analytics.js")),
  # useShinyjs(),
  # extendShinyjs(text = jsCode), #, functions = c("winprint")),
  
  useShinyalert(),
  
  # chooseSliderSkin(skin = "Shiny", color = "#e46c0a"),
  
  # chooseSliderSkin(skin = "Round", color = "#e46c0a"),
  
  # sandstone, spacelab Shiny themes
  navbarPage(theme = shinytheme("united"),
             id = "gMCStool", 
             title = "gMCStool",
             position = c("static-top"),
             # title = div(img(src = 'logo2.jpg',height = 28,width = 28), "gMCStool"),
             selected = "Overview", # selected = "1. Select Cell Lines", # (testing) !!!!!!!!!
             
             
             # 0: Overview
             source("ui-tab-overview.R",local=TRUE)$value, # source executes and loads the .R file
             
             # 1: HELP
             source("ui-tab-download_gMCS_lists.R",local=TRUE)$value,
             
             # 2: UPLOAD THE DATA
             source("ui-tab-upload_RNAseq_data.R",local=TRUE)$value,
             
             # 3: FIND ESSENTIAL GENES
             source("ui-tab-predict_essential_gene.R",local=TRUE)$value,
             
             # 4: VISUALIZE CASE-BY-CASE
             source("ui-tab-visualization.R",local=TRUE)$value,
             
             # 5: VISUALIZE CASE-BY-CASE
             source("ui-tab-DepMap_analysis.R",local=TRUE)$value,
             
             # 6: HELP
             source("ui-tab-help.R",local=TRUE)$value,
             
             # 7: ABOUT
             source("ui-tab-about.R",local=TRUE)$value
  ),
  
  # Footer
  source("ui-footer.R",local=TRUE)$value
)


########################################################################################################
#  SERVER
########################################################################################################

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # browser()
  
  
  # hide all the divs for hidden variables
  hideElement("mp2")
  hideElement("mp2_1");hideElement("mp2_2");hideElement("mp2_3");hideElement("mp2_4")
  hideElement("mp2_loading")
  hideElement("mp3")
  hideElement("mp3_loading")
  hideElement("mp4")
  hideElement("mp4_heatmap")
  hideElement("mp4_heatmap_1"); hideElement("mp4_heatmap_2"); hideElement("mp4_heatmap_3"); hideElement("mp4_loading")
  hideElement("mp5"); hideElement("mp5_dotplot_1"); hideElement("mp5_dotplot_2"); hideElement("mp5_dotplot_3")
  showElement("mp5")
  
  hideElement("mp_gMCS_GMCS_LIST_max_length") # to be included later on with the paper of Danel
  
  # generate hidden variables
  flags <- reactiveValues(flag_show_summary_data = F,
                          flag_show_table_essentiality = F,
                          flag_show_table_gmcs_essentiality = F,
                          flag_show_heatmap_gmcs = F,
                          flag_show_table_gmcs_correlation = F,
                          flag_show_dotplot_gmcs = F,
                          flag_show_images_gmcs_database = F)
  
  
  # generate variables ####
  values <- reactiveValues(I_GMCS_LIST = "EssentialTasks_CultureMedium",
                           gene.exp = NULL,
                           sample.class = as.factor("---"),
                           sample.cohort = as.factor("---"),
                           ResultsEssentiality = NULL,
                           filtersGEA = lapply(1:num_max_classes,function(x){c(0,100)}),
                           list.gMCS.essential.filtered = NULL,
                           list.corr.essential.filtered = NULL,
                           sample.class.target = rep(F,num_max_classes),
                           plot.obj = list(ggplot() + theme_void(), ggplot() + theme_void()),
                           plot.obj.DepMap = ggplot() + theme_void(),
                           O_table_essential_gmcs_rows_selected_previous_selection = 0,
                           O_table_essential_gmcs_bis_rows_selected_previous_selection = 0,
                           O_txtout_I_GMCS_LIST  = NULL,
                           O_txtout_I_TH_METHOD = NULL,
                           plot.5.piechart.gmcs = ggplot() + theme_void(),
                           plot.5.barchart.merged.gmcss = ggplot() + theme_void(),
                           plot.5.barchart.gmcss = ggplot() + theme_void(),
                           idx.messages = c(),
                           DepMap.info.filtered = DepMap.info.all[c("dictionary.CCLE")],
                           selected_tasks_previous = metTasks$EssentialTasks_CultureMedium$DESCRIPTION,
                           selected_gmcs_length_previous = max(gMCS.info.raw$EssentialTasks_CultureMedium$gMCSs.ENSEMBL.length)
                           
  )
  
  # change the gMCS info to be able to change in the code the custom
  gMCS.info <- reactiveValues(selected = gMCS.info.raw$EssentialTasks_CultureMedium)
  
  metTasks <- reactiveValues(EssentialTasks_CultureMedium = metTasks$EssentialTasks_CultureMedium,
                             EssentialTasks_FullMedium = metTasks$EssentialTasks_FullMedium, 
                             Only_CultureMedium = metTasks$Only_CultureMedium, 
                             Only_FullMedium = metTasks$Only_FullMedium, 
                             Custom_CultureMedium = metTasks$Custom_CultureMedium, 
                             Custom_FullMedium = metTasks$Custom_FullMedium)
  
  
  options(shiny.maxRequestSize = nBytesRAM) # 2 Gb of RAM
  
  # close all notifications with panel change
  observeEvent(input$gMCStool, {
    # print(paste("input$gMCStool", input$gMCStool))
    lapply(values$idx.messages, removeNotification)
    values$idx.messages <- c()
  })
  
  ## Buttons to change the page ####
  # 1) Change from overview to Data Upload
  observeEvent(c(input$I_mp0_next),{
    updateTabsetPanel(session, "gMCStool", selected = "2. Upload RNA-seq data")
  },ignoreInit = TRUE)
  observeEvent(c(input$I_mp0_next_to_download),{
    updateTabsetPanel(session, "gMCStool", selected = "1. gMCS database")
  },ignoreInit = TRUE)
  # 1) Change from Data Upload to Gene Essentially Analysis
  observeEvent(c(input$I_mp2_next),{
    updateTabsetPanel(session, "gMCStool", selected = "3. Predict Essential Genes")
  },ignoreInit = TRUE)
  # 1) Change from Gene Essentially Analysis to Case by Case
  observeEvent(c(input$I_mp3_next),{
    updateTabsetPanel(session, "gMCStool", selected = "4. Visualization")
  },ignoreInit = TRUE)
  
  
  # auxiliar function to hide panels at once
  hide_all_elements <- function(avoid = c()){
    tohide = c("mp2", "mp2_1", "mp2_2", "mp2_3", "mp2_4", 
               "mp3","mp3_loading", 
               "mp4", "mp4_heatmap_1", "mp4_heatmap_2", "mp4_heatmap_3",  
               "mp5_dotplot_1", "mp5_dotplot_2", "mp5_dotplot_3")
    
    for(k in setdiff(tohide, avoid)){hideElement(k)}
    
    flags$flag_show_summary_data <- F
    flags$flag_show_table_essentiality <- F
    flags$flag_show_table_gmcs_essentiality <- F
    flags$flag_show_table_gmcs_correlation <- F
    flags$flag_show_heatmap_gmcs <- F
    flags$flag_show_dotplot_gmcs <- F
    
  }
  
  # browser()
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # 1: gMCS database ####
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  
  # source("server-tab-download_gMCS_lists.R")
  
  # initialize tree
  output$I_gMCS_GMCS_LIST_TASKS_tree <- renderTree({
    lapply(split(metTasks[[1]]$DESCRIPTION, factor(metTasks[[1]]$group, levels = unique(metTasks[[1]]$group))),
           function(z){zz <- as.list(z); names(zz) <- z; return(structure(zz, stselected=TRUE))})    })
  
  
  # multiple selection of database
  observeEvent(c(input$I_GMCS_2_EssentialTask_or_Only, input$I_GMCS_2_CultureMedia_or_Full), {
    values$I_GMCS_LIST <- paste0(input$I_GMCS_2_EssentialTask_or_Only, '_', input$I_GMCS_2_CultureMedia_or_Full)
    
    gMCS.info$selected <- gMCS.info.raw[[values$I_GMCS_LIST]]
    
    updateRadioButtons(session = session, inputId = "I_GMCS_1_EssentialTask_or_Only", selected = input$I_GMCS_2_EssentialTask_or_Only)
    updateRadioButtons(session = session, inputId = "I_GMCS_1_CultureMedia_or_Full", selected = input$I_GMCS_2_CultureMedia_or_Full)
    
    values$O_txtout_I_GMCS_LIST <- gMCS.info$selected$fullname
    
    if (input$I_GMCS_2_EssentialTask_or_Only != "Custom") {
      updateSelectInput(session = session, inputId = "I_DepMap_GMCS_LIST", selected = values$I_GMCS_LIST)
    }
    
    # update slider of maximum length
    values$selected_tasks_previous <- max(gMCS.info.raw[[values$I_GMCS_LIST]]$gMCSs.ENSEMBL.length)
    updateSliderInput(session = session, inputId = "I_gMCS_GMCS_LIST_max_length",
                      value = max(gMCS.info.raw[[values$I_GMCS_LIST]]$gMCSs.ENSEMBL.length),
                      max = max(gMCS.info.raw[[values$I_GMCS_LIST]]$gMCSs.ENSEMBL.length), min = 0)
    
  })
  
  # Select all and unselect all ####
  observeEvent(input$I_gMCS_GMCS_LIST_TASKS_all,{
    # select all the tasks
    # use all the tasks as template
    idx <- paste0("EssentialTasks", '_', input$I_GMCS_2_CultureMedia_or_Full)
    
    output$I_gMCS_GMCS_LIST_TASKS_tree <- renderTree({
      lapply(split(metTasks[[idx]]$DESCRIPTION, factor(metTasks[[idx]]$group, levels = unique(metTasks[[idx]]$group))),
             function(z){zz <- as.list(z); names(zz) <- z; return(structure(zz, stselected=TRUE))})    })
  })
  
  observeEvent(input$I_gMCS_GMCS_LIST_TASKS_none,{
    # unselect all the metabolic tasks and show an empty table
    idx <- paste0("EssentialTasks", '_', input$I_GMCS_2_CultureMedia_or_Full)
    output$I_gMCS_GMCS_LIST_TASKS_tree <- renderTree({
      lapply(split(metTasks[[idx]]$DESCRIPTION, factor(metTasks[[idx]]$group, levels = unique(metTasks[[idx]]$group))),
             function(z){zz <- as.list(z); names(zz) <- z; return(structure(zz, stselected=FALSE))})    })
    
    # show empty table and empty figures
    output$O_gMCS_metabolic_tasks <- DT::renderDataTable({NULL} )
    values$plot.5.piechart.gmcs <- ggplot() + theme_void()
    values$plot.5.barchart.merged.gmcss <- ggplot() + theme_void()
    values$plot.5.barchart.gmcss <- ggplot() + theme_void()
    output$O_gMCS_piechart_gmcs <- renderPlot({values$plot.5.piechart.gmcs})
    output$O_gMCS_barchart_merged_gmcs <- renderPlot({values$plot.5.barchart.merged.gmcss})
    output$O_gMCS_barchart_gmcs <- renderPlot({values$plot.5.barchart.gmcss})
  })
  
  
  # Update the table
  observeEvent(values$I_GMCS_LIST, {
    #   # browser()
    idx <- ifelse(grepl("CultureMedium", values$I_GMCS_LIST), "EssentialTasks_CultureMedium", "EssentialTasks_FullMedium")
    if (grepl("Only", values$I_GMCS_LIST)) {
      idx2 <- grep("biomass", metTasks[[idx]]$DESCRIPTION)
      values$selected_tasks_previous <- metTasks[[idx]]$DESCRIPTION[idx2]
      
      output$I_gMCS_GMCS_LIST_TASKS_tree <- renderTree({
        c(lapply(split(metTasks[[idx]]$DESCRIPTION[-idx2], factor(metTasks[[idx]]$group[-idx2], levels = unique(metTasks[[idx]]$group[-idx2]))),
                 function(z){zz <- as.list(z); names(zz) <- z; return(structure(zz, stselected=FALSE))}),
          lapply(split(metTasks[[idx]]$DESCRIPTION[idx2], factor(metTasks[[idx]]$group[idx2], levels = unique(metTasks[[idx]]$group[idx2]))),
                 function(z){zz <- as.list(z); names(zz) <- z; return(structure(zz, stselected=T))}))
      })
    } else if ((grepl("EssentialTasks", values$I_GMCS_LIST))) {
      output$I_gMCS_GMCS_LIST_TASKS_tree <- renderTree({
        lapply(split(metTasks[[idx]]$DESCRIPTION, factor(metTasks[[idx]]$group, levels = unique(metTasks[[idx]]$group))),
               function(z){zz <- as.list(z); names(zz) <- z; return(structure(zz, stselected=TRUE))})    })
      values$selected_tasks_previous <- metTasks[[idx]]$DESCRIPTION
    }
  })
  
  ##  generate Custom task gMCS list ####
  observeEvent(c(input$I_gMCS_GMCS_LIST_TASKS_tree), {
    # store all
    all_tasks <- ifelse(grepl("CultureMedium", values$I_GMCS_LIST), "EssentialTasks_CultureMedium", "EssentialTasks_FullMedium")
    all_tasks <- metTasks[[all_tasks]]$DESCRIPTION
    
    # selected_tasks <- input$I_gMCS_GMCS_LIST_TASKS
    tree <- input$I_gMCS_GMCS_LIST_TASKS_tree
    req(tree)
    selected_tasks <- unlist(get_selected(tree))[grepl("^\\[", unlist(get_selected(tree)))]
    
    if (length(values$selected_tasks_previous) != length(selected_tasks) |
        any(values$selected_tasks_previous!=selected_tasks)){
      
      if (all(all_tasks %in% selected_tasks)){
        updateRadioButtons(session = session, inputId = "I_GMCS_1_EssentialTask_or_Only", selected = "EssentialTasks")
        updateRadioButtons(session = session, inputId = "I_GMCS_2_EssentialTask_or_Only", selected = "EssentialTasks")
      } else if (length(selected_tasks)==1 & sum(grepl("biomass production", selected_tasks))>0 ) {
        updateRadioButtons(session = session, inputId = "I_GMCS_1_EssentialTask_or_Only", selected = "Only")
        updateRadioButtons(session = session, inputId = "I_GMCS_2_EssentialTask_or_Only", selected = "Only")
      } else {
        updateRadioButtons(session = session, inputId = "I_GMCS_1_EssentialTask_or_Only", selected = "Custom")
        updateRadioButtons(session = session, inputId = "I_GMCS_2_EssentialTask_or_Only", selected = "Custom")
        
        values$I_GMCS_LIST <- paste0("Custom", '_', input$I_GMCS_2_CultureMedia_or_Full)
        
        origin_I_GMCS_LIST <- paste0("EssentialTasks", '_', input$I_GMCS_2_CultureMedia_or_Full)
        
        # reduce the number of tasks
        idx1 <- metTasks[[origin_I_GMCS_LIST]]$DESCRIPTION %in% selected_tasks
        metTasks[[values$I_GMCS_LIST]] <- metTasks[[origin_I_GMCS_LIST]][idx1,]
        
        # index of the gMCS that are included
        # browser()
        idx2 <- gMCS.info.raw[[origin_I_GMCS_LIST]]$table.gMCSs %>% as.data.frame() %>% 
          filter(task %in% which(idx1)) %>% 
          dplyr::select("idx") %>% as.matrix() %>% as.numeric() %>% unique() %>% sort()
        # idx2 <- sort(unique(gMCS.info[[origin_I_GMCS_LIST]]$table.gMCSs$idx[gMCS.info[[origin_I_GMCS_LIST]]$table.gMCSs$task %in% which(idx)]))
        
        # update gMCS info
        gMCS.info$selected <- as.list(gMCS.info.raw[[origin_I_GMCS_LIST]])
        gMCS.info$selected$gMCSs.ENSEMBL.txt <- gMCS.info.raw[[origin_I_GMCS_LIST]]$gMCSs.ENSEMBL.txt[idx2]
        gMCS.info$selected$gMCSs.ENSEMBL.length <- gMCS.info.raw[[origin_I_GMCS_LIST]]$gMCSs.ENSEMBL.length[idx2]
        # gMCS.info$selected$gMCSs.ENSEMBL.list <- gMCS.info.raw[[origin_I_GMCS_LIST]]$gMCSs.ENSEMBL.list[idx2]
        gMCS.info$selected$gMCSs.ENSEMBL.mat <- gMCS.info.raw[[origin_I_GMCS_LIST]]$gMCSs.ENSEMBL.mat[idx2,]
        gMCS.info$selected$gMCSs.ENSEMBL.txt.SYMBOL <- gMCS.info.raw[[origin_I_GMCS_LIST]]$gMCSs.ENSEMBL.txt.SYMBOL[idx2]
        gMCS.info$selected$gMCSs.ENSEMBL <- gMCS.info.raw[[origin_I_GMCS_LIST]]$gMCSs.ENSEMBL[idx1]
        
        # change a little bit the gMCS table for the new index
        gMCS.info$selected$table.gMCSs <- gMCS.info.raw[[origin_I_GMCS_LIST]]$table.gMCSs %>% filter(task %in% which(idx1))
        gMCS.info$selected$table.gMCSs$idx <- match(gMCS.info$selected$table.gMCSs$idx, idx2)
      }
      # save the results to avoid recalculating
      values$selected_tasks_previous <- selected_tasks
    }
    
    # update the slider for maximum length
    values$selected_tasks_previous <- max(gMCS.info.raw[[values$I_GMCS_LIST]]$gMCSs.ENSEMBL.length)
    updateSliderInput(session = session, inputId = "I_gMCS_GMCS_LIST_max_length",
                      value = max(gMCS.info.raw[[values$I_GMCS_LIST]]$gMCSs.ENSEMBL.length),
                      max = max(gMCS.info.raw[[values$I_GMCS_LIST]]$gMCSs.ENSEMBL.length), min = 0)
    
    
  })
  
  
  ##  generate reduce the number of gMCS according to the length list ####
  observeEvent(c(input$I_gMCS_GMCS_LIST_max_length), {
    
    selected_gmcs_length <- input$I_gMCS_GMCS_LIST_max_length
    # browser()
    
    if ( values$selected_gmcs_length_previous != selected_gmcs_length){
      
      # if (selected_gmcs_length <=10) {browser()}
      # browser()
      # index of the gMCS that are included by length
      idx2 <- gMCS.info$selected$gMCSs.ENSEMBL.length <= selected_gmcs_length
      
      if (!all(idx2)){
        # update gMCS info
        print(sapply(gMCS.info$selected$gMCSs.ENSEMBL,dim))
        gMCS.info$selected$gMCSs.ENSEMBL <- lapply(gMCS.info$selected$gMCSs.ENSEMBL, function(x){
          # x <- gMCS.info$selected$gMCSs.ENSEMBL[['57']]
          z <- x[apply(x!="",1,sum)<=selected_gmcs_length,]
          if(!is.matrix(z) ){ z <- matrix(z, nrow = 1)}
          z <- z[,apply(x!="",2,sum)>0]
          if(!is.matrix(z) ){ z <- matrix(z, nrow = 1)}
          if(prod(dim(z))==0){z <- matrix(NA, nrow = 0, ncol = 0)}
          if(is.null(z)){z <- matrix(NA, nrow = 0, ncol = 0)}
          return(z)})
        print(sapply(gMCS.info$selected$gMCSs.ENSEMBL,dim))
        gMCS.info$selected$gMCSs.ENSEMBL.txt <- gMCS.info$selected$gMCSs.ENSEMBL.txt[idx2]
        gMCS.info$selected$gMCSs.ENSEMBL.length <- gMCS.info$selected$gMCSs.ENSEMBL.length[idx2]
        # gMCS.info$selected$gMCSs.ENSEMBL.list <- gMCS.info$selected$gMCSs.ENSEMBL.list[idx2]
        gMCS.info$selected$gMCSs.ENSEMBL.mat <- gMCS.info$selected$gMCSs.ENSEMBL.mat[idx2,]
        gMCS.info$selected$gMCSs.ENSEMBL.txt.SYMBOL <- gMCS.info$selected$gMCSs.ENSEMBL.txt.SYMBOL[idx2]
        
        # change a little bit the gMCS table for the new index
        gMCS.info$selected$table.gMCSs$idx <- match(gMCS.info$selected$table.gMCSs$idx, which(idx2))
        
      }
      # save the results to avoid recalculating
      values$selected_gmcs_length_previous <- selected_gmcs_length
    }
    
  })
  
  
  # generate table with all the selected tasks ####
  Selected_metabolic_tasks <- reactive({
    
    tree <- input$I_gMCS_GMCS_LIST_TASKS_tree
    req(tree)
    selected_tasks <- unlist(get_selected(tree))[grepl("^\\[", unlist(get_selected(tree)))]
    
    Selected_metabolic_tasks <- data.frame("task number" = metTasks[[values$I_GMCS_LIST]][,1],
                                           "task name" = metTasks[[values$I_GMCS_LIST]][,2],
                                           "task group" = metTasks[[values$I_GMCS_LIST]][,3],
                                           "num gMCSs" = unlist(lapply(gMCS.info$selected$gMCSs.ENSEMBL,nrow)))
    rownames(Selected_metabolic_tasks) <- 1:nrow(Selected_metabolic_tasks)
    
    Selected_metabolic_tasks <- Selected_metabolic_tasks[ Selected_metabolic_tasks[,"task.name"] %in% selected_tasks,]
    
    
    
    COLORS <- scales::hue_pal()(length(unique(Selected_metabolic_tasks[,"task.group"])))
    COLORS <- colorspace::adjust_transparency(COLORS, alpha = 0.1)
    names(COLORS) <- unique(Selected_metabolic_tasks[,"task.group"])
    
    COLORS <- COLORS[!is.na(names(COLORS))]
    
    Selected_metabolic_tasks <- datatable(Selected_metabolic_tasks,
                                          filter = "none", selection = "none", rownames = F,
                                          options = list(pageLength = -1,
                                                         info = FALSE,
                                                         lengthMenu = list(c(-1, 10, 20, 50, 100), c("All", "10", "20", "50", "100")))
    ) %>% formatStyle(columns = "task.group",
                      backgroundColor = styleEqual(names(COLORS), COLORS))
  })
  
  
  observeEvent(input$I_show_summary_images_gMCS_database,{
    # print(input$I_show_summary_images_gMCS_database)
    flags$flag_show_images_gmcs_database <- !flags$flag_show_images_gmcs_database
  })
  
  
  output$O_gMCS_piechart_gmcs <- renderPlot({values$plot.5.piechart.gmcs})
  output$O_gMCS_barchart_merged_gmcs <- renderPlot({values$plot.5.barchart.merged.gmcss})
  output$O_gMCS_barchart_gmcs <- renderPlot({values$plot.5.barchart.gmcss})
  
  # show the selected gMCS for certain gene
  observeEvent(c(input$I_gMCS_GMCS_LIST_TASKS_tree,
                 flags$flag_show_images_gmcs_database), {
                   
                   print(paste("flags$flag_show_images_gmcs_database", flags$flag_show_images_gmcs_database))
                   hideElement("mp_gMCS_1"); hideElement("mp_gMCS_2"); showElement("mp_gMCS_3"); showElement("mp_gMCS_download")
                   
                   tree <- input$I_gMCS_GMCS_LIST_TASKS_tree
                   req(tree)
                   selected_tasks <- unlist(get_selected(tree))[grepl("^\\[", unlist(get_selected(tree)))]
                   
                   Selected.metabolic.tasks <- metTasks[[values$I_GMCS_LIST]]
                   Selected.metabolic.tasks[,"num gMCSs"] <- unlist(lapply(gMCS.info$selected$gMCSs.ENSEMBL,nrow))
                   Selected.metabolic.tasks <- Selected.metabolic.tasks[Selected.metabolic.tasks[,2] %in% selected_tasks,]
                   
                   if (flags$flag_show_images_gmcs_database){
                     if (grepl("Only", values$I_GMCS_LIST, ignore.case = T)){
                       hideElement("mp_gMCS_1")
                       showElement("mp_gMCS_2")
                       values$plot.5.piechart.gmcs <- ggplot() + theme_void()
                       values$plot.5.barchart.merged.gmcss <- ggplot() + theme_void()
                       
                       # browser()
                       
                       sdf <- reshape2::melt(lapply(gMCS.info$selected$gMCSs.ENSEMBL, function(x){table(apply(x!="",1,sum))}))
                       colnames(sdf) <- c("num_genes", "num_gMCSs", "ID")
                       sdf$ID[sdf$ID==58] <- 57
                       
                       
                       sdf <- merge(Selected.metabolic.tasks, sdf, all = T)
                       sdf[is.na(sdf)] <- 0
                       # sdf2 <- sdf[input$I_gMCS_GMCS_LIST_TASKS %in% sdf[,2],]
                       
                       sdf2 <- lapply(selected_tasks, function(x){sdf[sdf$DESCRIPTION==x,]})
                       values$plot.5.barchart.gmcss <- lapply(sdf2,  function(x){ ggbarplot(x, x = "num_genes", y = "num_gMCSs",
                                                                                            facet.by = "ID", fill = "num_genes") + xlab("") + ylab("") +
                           facet_wrap(facets = vars(ID), scales = "free", nrow = 1, ncol = 1) +
                           theme(legend.position = "none")})
                     } else {
                       showElement("mp_gMCS_1")
                       hideElement("mp_gMCS_2")
                       values$plot.5.barchart.gmcss <- ggplot() + theme_void()
                       tab5 <- aggregate(Selected.metabolic.tasks$`num gMCSs`, by = list(DESCRIPTION = Selected.metabolic.tasks$group), sum)
                       tab5$DESCRIPTION <- factor(tab5$DESCRIPTION, levels = unique(Selected.metabolic.tasks$group))
                       tab5 <- tab5[order(tab5$DESCRIPTION,decreasing = F),]
                       tab5$DESCRIPTION_2 <- paste0(tab5$DESCRIPTION, ' (' ,tab5$x, ')')
                       tab5$DESCRIPTION_2 <- factor(tab5$DESCRIPTION_2, levels = tab5$DESCRIPTION_2)
                       
                       tab5$aux1 <- 1
                       
                       kkk <- c(0, cumsum(tab5$x))
                       idx <- sum(tab5$x) - unlist(lapply(1:dim(tab5)[1], function(i){median(kkk[c(i,i+1)])}))
                       
                       tab5$aux2 <- idx
                       tab5$aux3 <- 1.6
                       
                       values$plot.5.piechart.gmcs <- ggplot(tab5, aes(x = aux1, y = x, label = x, fill = DESCRIPTION_2, width = aux1)) +
                         geom_bar(stat = "identity", position = "stack", col = "white")  +
                         theme_bw() + #scale_fill_manual(values = colorRampPalette(palette)(length(levels(df$Var1)))) +
                         ggrepel::geom_label_repel(aes(y = aux2, x = aux3, fill = NULL),
                                                   nudge_x = 0,
                                                   nudge_y = 0,
                                                   show.legend = F, direction = "both") +
                         coord_polar("y", start=0, direction = -1) +
                         # ggtitle(paste0("Distribution of MM specific fusion genes per biotype, read throughs\nonlySpanningJunction ", date())) +
                         theme(panel.background = element_rect(color = "white", fill = "white"),
                               panel.border = element_rect(fill=NA, colour = "white"),
                               panel.grid = element_blank(),
                               axis.ticks = element_blank(),
                               # legend.margin = unit(10, 'mm'),
                               axis.text = element_blank(),
                               axis.title = element_blank(),
                               plot.background = element_blank(),
                               legend.title = element_text(colour = "white"),
                               plot.margin = unit(c(0,0,0,0), "mm"),
                               legend.position = "right")
                       
                       idx <- metTasks[[values$I_GMCS_LIST]]$DESCRIPTION %in% selected_tasks
                       if (all(idx)){
                         sdf <- table(gMCS.info$selected$gMCSs.ENSEMBL.length)
                       } else {
                         idx <- sort(unique(gMCS.info$selected$EssentialTasks_CultureMedium$table.gMCSs$idx[gMCS.info$selected$EssentialTasks_CultureMedium$table.gMCSs$task %in% which(idx)]))
                         sdf <- table(gMCS.info$selected$gMCSs.ENSEMBL.length[idx])
                       }
                       
                       sdf <- reshape2::melt(sdf)
                       colnames(sdf) <- c("gmcs.length", "Freq")
                       
                       
                       sdf2 <- sdf[sdf$gmcs.length<20, ]
                       sdf2$gmcs.length <- as.character(sdf2$gmcs.length)
                       sdf2 <- as.data.frame(rbind(sdf2,0))
                       sdf2[dim(sdf2)[1],1] <- ">20"
                       sdf2[dim(sdf2)[1],2] <- sum(sdf[sdf$gmcs.length>=20, 2])
                       sdf2$gmcs.length <- factor(sdf2$gmcs.length,levels = sdf2$gmcs.length)
                       
                       
                       values$plot.5.barchart.merged.gmcss <- ggplot(sdf2, aes(x = gmcs.length,y = Freq)) +
                         geom_bar(stat="identity", fill = rgb(235,110,4,255,maxColorValue = 255),color="black") +
                         coord_flip() +  theme_classic()+ xlab("Number of genes in the gMCSs") + ylab("Number of gMCSs") +
                         # scale_y_continuous(trans = "log10")+
                         geom_text(stat='identity', aes(label=Freq),size = 6, vjust=0.5, hjust = -0.5) +
                         theme(text = element_text(size = 16, face = "bold"),
                               axis.text = element_text(size = 16, face = "bold"),
                               axis.text.x = element_text(angle = 75, vjust=1, hjust = 1),
                               panel.background = element_blank() # bg of the panel
                               , plot.background = element_blank() # bg of the plot
                               , panel.grid.major = element_blank() # get rid of major grid
                               , panel.grid.minor = element_blank() # get rid of minor grid
                               , legend.background = element_blank() # get rid of legend bg
                               , legend.box.background = element_blank() # get rid of legend panel bg
                         )
                       if (any(grepl("growth", selected_tasks, ignore.case = T))) {
                         values$plot.5.barchart.merged.gmcss <- values$plot.5.barchart.merged.gmcss +
                           scale_y_continuous(limits = c(-100,50e3), expand = c(0, 0),breaks = seq(0,40e3,5e3), minor_breaks = seq(0,25e3,5e3))
                       }
                       
                     }
                   } else {
                     hideElement("mp_gMCS_1")
                     hideElement("mp_gMCS_2")
                     values$plot.5.piechart.gmcs <- ggplot() + theme_void()
                     values$plot.5.barchart.merged.gmcss <- ggplot() + theme_void()
                     values$plot.5.barchart.gmcss <- ggplot() + theme_void()
                   }
                   
                   output$O_gMCS_piechart_gmcs <- renderPlot({values$plot.5.piechart.gmcs})
                   output$O_gMCS_barchart_merged_gmcs <- renderPlot({values$plot.5.barchart.merged.gmcss})
                   output$O_gMCS_barchart_gmcs <- renderPlot({values$plot.5.barchart.gmcss})
                   output$O_gMCS_metabolic_tasks <- DT::renderDataTable({Selected_metabolic_tasks()} )
                 })
  
  output$O_gMCS_database_download <- downloadHandler(filename = paste0("gMCStool_selected_gMCSs_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".Rdata"), 
                                                     content = function(fname){gMCS.info <- gMCS.info$selected; save(gMCS.info, file = fname)})
  
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # 2: load the RNA-seq data ####
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  {
    # select menu to show
    observeEvent(input$I_DataInput_mode, {
      hide_all_elements()
      print(input$I_DataInput_mode)
      if (input$I_DataInput_mode=="text"){ showElement("mp2_1")
      } else if (input$I_DataInput_mode=="tximport"){ showElement("mp2_2")
      } else if (input$I_DataInput_mode=="Rdata"){ showElement("mp2_3")}
    })
    
    # select to show examples
    observeEvent(input$I_Data_input_examples, {
      hide_all_elements()
      print("I_Data_input_examples")
      showElement("mp2_4")
      updateRadioButtons(session, inputId = "I_DataInput_mode", selected = character(0))
    })
    
    
    
    # load the results
    observeEvent(input$Data_input_rdata,{
      if (is.null(input$Data_input_rdata) == FALSE){
        hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
        showElement("mp2_loading")
        flags$flag_show_summary_data <- F
        flags$flag_show_table_essentiality <- F
        flags$flag_show_heatmap_gmcs <- F
        flags$flag_show_dotplot_gmcs <- F
        
        print(input$Data_input_rdata)
        # browser()
        load(input$Data_input_rdata$datapath)
        
        values$gene.exp  <- gene.exp
        values$sample.class <- sample.class
        values$sample.cohort <- sample.cohort
        
        flags$flag_show_summary_data <- T
        hideElement("mp2_loading")
        showElement("mp2")
      }
    })
    
    # load precomputed results from example dataset
    observeEvent(input$I_Data_input_example_TPM,{
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp2_loading")
      
      load("ExampleData/ExampleData_HumanGEM_BcellMM_TPM.Rdata")
      
      values$gene.exp  <- gene.exp
      values$sample.class <- sample.class
      values$sample.cohort <- as.factor(rep("BcellMM", length(sample.class)))
      
      flags$flag_show_summary_data <- T
      hideElement("mp2_loading")
      showElement("mp2")
    })
    observeEvent(input$I_Data_input_example_log2TPM,{
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp2_loading")
      
      load("ExampleData/ExampleData_HumanGEM_BcellMM_log2TPM.Rdata")
      
      values$gene.exp  <- gene.exp
      values$sample.class <- sample.class
      values$sample.cohort <- as.factor(rep("BcellMM", length(sample.class)))
      
      flags$flag_show_summary_data <- T
      hideElement("mp2_loading")
      showElement("mp2")
    })
    
    observeEvent(input$I_Data_input_example_DepMap,{
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp2_loading")
      
      load("ExampleData/CompleteResults_HumanGEM_DepMap_gmcsTH5_TPM.Rdata")
      
      values$gene.exp  <- gene.exp
      values$sample.class <- sample.class
      values$sample.cohort <- sample.cohort
      
      flags$flag_show_summary_data <- T
      hideElement("mp2_loading")
      showElement("mp2")
    })
    
    # load the previously saved data
    observeEvent(input$Data_input_rdata,{
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp2_loading")
      
      print(input$Data_input_rdata$datapath)
      load(input$Data_input_rdata$datapath)
      
      values$gene.exp  <- gene.exp
      values$sample.class <- sample.class
      values$sample.cohort <- sample.cohort
      
      flags$flag_show_summary_data <- T
      hideElement("mp2_loading")
      showElement("mp2")
    })
    
    # read the table with the sample classification
    observeEvent(input$Data_input_text_class,{
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp2_loading")
      flags$flag_show_summary_data <- F
      flags$flag_show_table_essentiality <- F
      flags$flag_show_heatmap_gmcs <- F
      flags$flag_show_dotplot_gmcs <- F
      
      print(input$Data_input_text_class$datapath)
      sample.class <- as.data.frame(data.table::fread(input$Data_input_text_class$datapath))
      
      if (!is.null(values$gene.exp)){
        if (nrow(sample.class)==ncol(values$gene.exp)) {
          if (ncol(sample.class)==3) {
            rownames(sample.class) <- sample.class[,1]
            sample.class <- sample.class[colnames(values$gene.exp),]
            sample.class <- sample.class[,-1]
          }
          values$sample.class <- factor(sample.class[,1], levels = unique(sample.class[,1]))
          values$sample.cohort <-  factor(sample.class[,2], levels = unique(sample.class[,2]))
        } else {
          values$idx.messages <- c(values$idx.messages, showNotification("It is not possible to use this file:\n
                       the number of rows and the number of samples in the gene expression does not match",
                                                                         duration = NULL, type = "error"))
        } 
        flags$flag_show_summary_data <- T
      }
      
      hideElement("mp2_loading")
      showElement("mp2")
    })
    
    # read the table with the sample classification
    observeEvent(input$Data_input_text_class_2,{
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp2_loading")
      flags$flag_show_summary_data <- F
      flags$flag_show_table_essentiality <- F
      flags$flag_show_heatmap_gmcs <- F
      flags$flag_show_dotplot_gmcs <- F
      
      print(input$Data_input_text_class_2$datapath)
      sample.class <- as.data.frame(data.table::fread(input$Data_input_text_class_2$datapath))
      
      
      if (nrow(sample.class)==ncol(values$gene.exp)) {
        if (ncol(sample.class)==3) {
          rownames(sample.class) <- sample.class[,1]
          sample.class <- sample.class[colnames(values$gene.exp),]
          sample.class <- sample.class[,-1]
        }
        values$sample.class <- factor(sample.class[,1], levels = unique(sample.class[,1]))
        values$sample.cohort <-  factor(sample.class[,2], levels = unique(sample.class[,2]))
      } else {
        values$idx.messages <- c(values$idx.messages, showNotification("It is not possible to use this file:\n
                       the number of rows and the number of samples in the gene expression does not match",
                                                                       duration = NULL, type = "error"))
      } 
      
      
      flags$flag_show_summary_data <- T
      hideElement("mp2_loading")
      showElement("mp2")
    })
    
    
    # sort the levels in sample.cohort and sample.class
    observeEvent(input$I_Data_input_alpha_sort_class,{
      values$sample.class <- factor(as.character(values$sample.class),
                                    levels = sort(unique(as.character(values$sample.class)), decreasing = input$I_Data_input_alpha_sort_class%%2==0))
    })
    observeEvent(input$I_Data_input_alpha_sort_cohort,{
      values$sample.cohort <- factor(as.character(values$sample.cohort),
                                     levels = sort(unique(as.character(values$sample.cohort)), decreasing = input$I_Data_input_alpha_sort_cohort%%2==0))
    })
    
    observeEvent(input$I_Data_input_num_sort_class,{
      values$sample.class <- factor(as.character(values$sample.class),
                                    levels = levels(values$sample.class)[order(table(values$sample.class), decreasing = input$I_Data_input_num_sort_class%%2==0)])
    })
    observeEvent(input$I_Data_input_num_sort_cohort,{
      values$sample.cohort <- factor(as.character(values$sample.cohort),
                                     levels = levels(values$sample.cohort)[order(table(values$sample.cohort), decreasing = input$I_Data_input_num_sort_cohort%%2==0)])
    })
    
    
    # relevel the tables if the sample.class or cohort has been reordered
    observeEvent(c(values$sample.class,
                   values$sample.cohort),{
                     
                     print("relevel the tables")
                     
                     if (!is.null(values$ResultsEssentiality)){
                       # browser()
                       list.gene.essential.aux <- values$ResultsEssentiality$num.essential.gene
                       idx <- which(colnames(list.gene.essential.aux) %in% levels(values$sample.class))
                       list.gene.essential.aux[,idx] <- list.gene.essential.aux[,levels(values$sample.class)]
                       colnames(list.gene.essential.aux)[idx] <- levels(values$sample.class)
                       values$ResultsEssentiality$num.essential.gene <- list.gene.essential.aux
                       
                       list.gene.essential.aux <- values$ResultsEssentiality$ratio.essential.gene
                       idx <- which(colnames(list.gene.essential.aux) %in% levels(values$sample.class))
                       list.gene.essential.aux[,idx] <- list.gene.essential.aux[,levels(values$sample.class)]
                       colnames(list.gene.essential.aux)[idx] <- levels(values$sample.class)
                       values$ResultsEssentiality$ratio.essential.gene <- list.gene.essential.aux
                       
                       list.gene.essential.aux <- values$ResultsEssentiality$list.gMCS.essential
                       idx <- which(colnames(list.gene.essential.aux) %in% levels(values$sample.class))
                       list.gene.essential.aux[,idx] <- list.gene.essential.aux[,levels(values$sample.class)]
                       colnames(list.gene.essential.aux)[idx] <- levels(values$sample.class)
                       values$ResultsEssentiality$list.gMCS.essential <- list.gene.essential.aux
                       
                     }
                   })
    
    
    
    # read the table with the sample classification
    observeEvent(input$Data_input_text_expr,{
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp2_loading")
      flags$flag_show_summary_data <- F
      flags$flag_show_table_essentiality <- F
      flags$flag_show_heatmap_gmcs <- F
      flags$flag_show_dotplot_gmcs <- F
      
      # browser()
      
      print(input$Data_input_text_expr$datapath)
      gene.exp <- as.data.frame(data.table::fread(input$Data_input_text_expr$datapath))
      rownames(gene.exp) <- gene.exp[,1]
      gene.exp <- gene.exp[,-1]
      
      # check the number of genes that are in the gMCS
      genes.gMCS <- unique(unlist(lapply(gMCS.info.raw[1:4], function(x){x$genes.gMCSs.ENSEMBL})))
      
      coverage <- mean(genes.gMCS %in% rownames(gene.exp))
      
      if (ncol(gene.exp)>num_max_samples){
        values$idx.messages <- c(values$idx.messages, showNotification(paste0("This dataset is too large for the online tool,\n",
                                                                              "please download the tool from GitHub"),
                                                                       duration = NULL, type = "error"))
        flags$flag_show_summary_data <- F
        hideElement("mp2")
      }
      
      if (coverage < 0.5) {
        values$idx.messages <- c(values$idx.messages, showNotification(paste0("Most of the genes of the gMCS are not in this dataset\n",
                                                                              "coverage = ", round(coverage*100,2), "%"),
                                                                       duration = NULL, type = "error"))
        flags$flag_show_summary_data <- F
        hideElement("mp2")
      } else {
        values$gene.exp <- gene.exp
        values$sample.class <- rep('--', ncol(gene.exp))
        values$sample.cohort <- rep('--', ncol(gene.exp))
        
        
        values$idx.messages <- c(values$idx.messages, showNotification(paste0("Most of the genes of the gMCS are not in this dataset\n",
                                                                              "coverage = ", round(coverage*100,2), "%"),
                                                                       duration = NULL, type = ifelse(coverage < 0.9, "warning", "message")))
        
        flags$flag_show_summary_data <- T
        showElement("mp2")
      }
      
      hideElement("mp2_loading")
    })
    
    
    
    
    # read the table with the sample classification
    observeEvent(input$Data_input_tximport,{
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp2_loading")
      flags$flag_show_summary_data <- F
      flags$flag_show_table_essentiality <- F
      flags$flag_show_heatmap_gmcs <- F
      flags$flag_show_dotplot_gmcs <- F
      
      print(input$Data_input_tximport$datapath)
      
      if (grepl("\\.rds$", input$Data_input_tximport$datapath, ignore.case = T)) {
        tximport <- readRDS(input$Data_input_tximport$datapath)
      } else if (grepl("\\.rdata$", input$Data_input_tximport$datapath, ignore.case = T)) {
        tximport <- new.env()
        load(input$Data_input_tximport$datapath, envir = tximport)
        tximport <- tximport[[1]]
      }  
      
      try({
        gene.exp <- tximport$abundance
        rownames(gene.exp) <- sub('\\..*', '', rownames(gene.exp))
        
        if (mean(grepl("ENSG", rownames(gene.exp))<0.5)) {
          showNotification(paste0("Genes are not in ENSEMBL (ENSG)"),
                           duration = NULL, type = "error")
        }
        
        # check the number of genes that are in the gMCS
        genes.gMCS <- unique(unlist(lapply(gMCS.info.raw[1:4], function(x){x$genes.gMCSs.ENSEMBL})))
        coverage <- mean(genes.gMCS %in% rownames(gene.exp))
        
        if (ncol(gene.exp)>num_max_samples){
          values$idx.messages <- c(values$idx.messages, showNotification(paste0("This dataset is too large for the online tool,\n",
                                                                                "please download the tool from GitHub"),
                                                                         duration = NULL, type = "error"))
          flags$flag_show_summary_data <- F
          hideElement("mp2")
        }
        
        if (coverage < 0.5) {
          values$idx.messages <- c(values$idx.messages, showNotification(paste0("Most of the genes of the gMCS are not in this dataset\n",
                                                                                "coverage = ", round(coverage*100,2), "%"), duration = NULL, type = "error"))
          flags$flag_show_summary_data <- F
          hideElement("mp2")
        } else {
          values$gene.exp <- gene.exp
          values$sample.class <- rep('--', ncol(gene.exp))
          values$sample.cohort <- rep('--', ncol(gene.exp))
          
          values$idx.messages <- c(values$idx.messages, showNotification(paste0("Most of the genes of the gMCS are not in this dataset\n",
                                                                                "coverage = ", round(coverage*100,2), "%"), duration = NULL, type = ifelse(coverage < 0.9, "warning", "message")))
          
          flags$flag_show_summary_data <- T
          showElement("mp2")
        }
      })
      
      hideElement("mp2_loading")
    })
    
    
    # buttons to download the data ####
    output$O_data_input_download_1 <- downloadHandler(filename = paste0("gMCStool_gene_expression_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".txt"), 
                                                      content = function(fname){gene.exp  <- values$gene.exp; fwrite(gene.exp, file = fname, quote = F, row.names = T, sep = ";")})
    
    output$O_data_input_download_2 <- downloadHandler(filename = paste0("gMCStool_sample_classification_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".txt"), 
                                                      content = function(fname){
                                                        zz <- data.frame(colnames(values$gene.exp), values$sample.class, values$sample.cohort)
                                                        colnames(zz) <- c("Sample.ID", "sample.class", "sample.cohort")
                                                        fwrite(zz, file = fname, quote = F, row.names = F, sep = ";")})
    
    output$O_data_input_download_rdata <- downloadHandler(filename = paste0("gMCStool_input_data_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".Rdata"), 
                                                          content = function(fname){gene.exp  <- values$gene.exp; sample.class <- values$sample.class; sample.cohort <- values$sample.cohort; save(gene.exp, sample.class, sample.cohort, file = fname)})
    
    
    
    
    # calculate which table to show ####
    SummaryClassCohortTable <- reactive({
      if (flags$flag_show_summary_data & !is.null(values$gene.exp)){
        zz <- data.frame(colnames(values$gene.exp), values$sample.class, values$sample.cohort)
        colnames(zz) <- c("Sample.ID", "sample.class", "sample.cohort")
        
        SummaryClassCohortTable <- datatable(zz,
                                             selection = "none", filter = "none",
                                             options = list(pageLength = -1, 
                                                            info = FALSE,
                                                            lengthMenu = list(c(-1, 10, 20, 50, 100), c("All", "10", "20", "50", "100"))),
                                             editable = T)
        
        if (length(levels(values$sample.class))>1){
          COLORS <- scales::hue_pal()(length(levels(values$sample.class)))
          names(COLORS) <- levels(values$sample.class)
          COLORS <- COLORS[!is.na(names(COLORS))]
          COLORS <- colorspace::adjust_transparency(COLORS, alpha = 0.1)
          names(COLORS) <- levels(values$sample.class)
          SummaryClassCohortTable <- SummaryClassCohortTable %>% formatStyle(columns = "sample.class",
                                                                             backgroundColor = styleEqual(names(COLORS), COLORS))
          # }
          # if (length(levels(values$sample.cohort))>1){
          COLORS2 <- viridis::viridis(length(levels(values$sample.cohort)))
          names(COLORS2) <- levels(values$sample.cohort)
          COLORS2 <- COLORS2[!is.na(names(COLORS2))]
          COLORS2 <- colorspace::adjust_transparency(COLORS2, alpha = 0.1)
          names(COLORS2) <- levels(values$sample.cohort)
          SummaryClassCohortTable <- SummaryClassCohortTable %>% formatStyle(columns = "sample.cohort",
                                                                             backgroundColor = styleEqual(names(COLORS2), COLORS2))
        }
        
        
        SummaryClassCohortTable <- SummaryClassCohortTable
        
      } else {
        hideElement("mp2")
        SummaryClassCohortTable <- NULL
      }
    })
    
    SummaryClassCohortTable_1 <- reactive({
      if (flags$flag_show_summary_data & !is.null(values$sample.class)){
        SummaryClassCohortTable_1 <- as.data.frame(t(as.matrix(table(values$sample.class))))
      } else { SummaryClassCohortTable_1 <- NULL }
    })
    SummaryClassCohortTable_2 <- reactive({
      if (flags$flag_show_summary_data & !is.null(values$sample.cohort)){
        SummaryClassCohortTable_2 <- as.data.frame(t(as.matrix(table(values$sample.cohort))))
      } else { SummaryClassCohortTable_2 <- NULL }
    })
    
    
    # summary tables and total table of samples
    observeEvent(flags$flag_show_summary_data, {
      print(flags$flag_show_summary_data)
      if (!flags$flag_show_summary_data){
        hideElement("mp2")
        hideElement("mp2_loading")
        output$O_table_sample_class_cohort <- NULL
        output$O_table_sample_class_cohort_1 <- NULL
        output$O_table_sample_class_cohort_2 <- NULL
      } else {
        showElement("mp2")
        output$O_table_sample_class_cohort <- DT::renderDataTable({SummaryClassCohortTable()},
                                                                  options = list(pageLength = -1, 
                                                                                 info = FALSE,
                                                                                 lengthMenu = list(c(-1, 10, 20, 50, 100), c("All", "10", "20", "50", "100"))))
        output$O_table_sample_class_cohort_1 <- renderTable({SummaryClassCohortTable_1()})
        output$O_table_sample_class_cohort_2 <- renderTable({SummaryClassCohortTable_2()})
      }
    })
    
    
    
    # edit the data of the table (works!) ####
    proxy = dataTableProxy("O_table_sample_class_cohort")
    
    observeEvent(input$O_table_sample_class_cohort_cell_edit, {
      info = input$O_table_sample_class_cohort_cell_edit
      str(info)
      i = info$row
      j = info$col
      k = info$value
      str(paste(i,j,k))
      
      isolate(
        if (j == 2) {
          sample.class <- as.character(values$sample.class)
          sample.class[i] <- DT::coerceValue(k, sample.class[i])
          values$sample.class <- factor(sample.class, levels = unique(sample.class))
        } else if (j == 3) {
          sample.cohort <- as.character(values$sample.cohort)
          sample.cohort[i] <- DT::coerceValue(k, sample.cohort[i])
          values$sample.cohort <- factor(sample.cohort, levels = unique(sample.cohort))
        } else {
          warning("You cannot change sample names.") # check to stop the user from editing only few columns
        }
      )
      flags$flag_show_summary_data <- T
      # replaceData(proxy, v$data, resetPaging = FALSE)  # replaces data displayed by the updated table
    })
    
  }
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # 3. FIND ESSENTIAL GENES ####
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  {
    # multiple selection of database
    observeEvent(c(input$I_GMCS_1_EssentialTask_or_Only, input$I_GMCS_1_CultureMedia_or_Full), {
      values$I_GMCS_LIST <- paste0(input$I_GMCS_1_EssentialTask_or_Only, '_', input$I_GMCS_1_CultureMedia_or_Full)
      updateRadioButtons(session = session, inputId = "I_GMCS_2_EssentialTask_or_Only", selected = input$I_GMCS_1_EssentialTask_or_Only)
      updateRadioButtons(session = session, inputId = "I_GMCS_2_CultureMedia_or_Full", selected = input$I_GMCS_1_CultureMedia_or_Full)
      
      values$O_txtout_I_GMCS_LIST <- gMCS.info$selected$fullname
      
      if (input$I_GMCS_1_EssentialTask_or_Only != "Custom") {
        updateSelectInput(session = session, inputId = "I_DepMap_GMCS_LIST", selected = values$I_GMCS_LIST)
      }
      
    })
    
    
    # menu to show in left column 
    observeEvent(input$I_TH_METHOD, {
      if (input$I_TH_METHOD=="localT2"){
        showElement("I_perc_localT2_show")
        hideElement("I_perc_gmcsTH_show")
        hideElement("I_perc_singleTH_show")
      } else if (input$I_TH_METHOD=="gmcsTH"){
        hideElement("I_perc_localT2_show")
        showElement("I_perc_gmcsTH_show")
        hideElement("I_perc_singleTH_show")
      } else if (input$I_TH_METHOD=="singleTH"){
        hideElement("I_perc_localT2_show")
        hideElement("I_perc_gmcsTH_show")
        showElement("I_perc_singleTH_show")
      } 
    })
    
    # sincronize text and slider 
    observeEvent(input$I_perc_gmcsTH_txt,{
      if(as.numeric(input$I_perc_gmcsTH_txt) != input$I_perc_gmcsTH){
        updateSliderInput(session = session,inputId = 'I_perc_gmcsTH',value = input$I_perc_gmcsTH_txt, step = 0.001)
      }
    })
    observeEvent(input$I_perc_gmcsTH,{
      if(as.numeric(input$I_perc_gmcsTH_txt) != input$I_perc_gmcsTH){
        updateTextInput(session = session,inputId = 'I_perc_gmcsTH_txt',value = input$I_perc_gmcsTH) 
        updateSliderInput(session = session,inputId = 'I_perc_gmcsTH', step = 0.01)
      }
    })
    
    
    observeEvent(input$I_action_predictEssentialGene, {
      print(paste("Calculate = ", input$I_action_predictEssentialGene))
      hideElement("mp3_ResultsEssentiality_input_RDATA")
      hideElement("mp3_ResultsEssentiality_input_example")
    })
    # run the calculation of essential genes ####
    observeEvent(input$I_action_predictEssentialGene,{
      # updateDataTable <<- 1
      # source("fun-CalculateEssentialGenes_gmcsTH.R")
      # browser()
      if (!is.null(values$gene.exp) & !is.null(values$sample.class) & !is.null(values$sample.cohort)){
        
        hideElement("mp3")
        showElement("mp3_loading")
        flags$flag_show_table_essentiality <- F
        
        if (input$I_TH_METHOD=="gmcsTH"){
          # browser()
          values$O_txtout_I_GMCS_LIST <- gMCS.info$selected$fullname
          values$O_txtout_I_TH_METHOD <- paste0("gmcsTH (", input$I_perc_gmcsTH*100, "%)")
          ResultsEssentiality <- CalculateEssentialGenes_gmcsTH(gene.exp = values$gene.exp,
                                                                gMCS.info = gMCS.info$selected,
                                                                sample.class = values$sample.class,
                                                                gmcsTH_perc = input$I_perc_gmcsTH, 
                                                                nWorkers = nWorkers)
        } else if (input$I_TH_METHOD=="localT2") {
          values$O_txtout_I_GMCS_LIST <- gMCS.info$selected$fullname
          values$O_txtout_I_TH_METHOD <- paste0("localT2 (", input$I_perc_localT2, ")")
          ResultsEssentiality <- CalculateEssentialGenes_localT2(gene.exp = values$gene.exp,
                                                                 gMCS.info = gMCS.info$selected,
                                                                 sample.class = values$sample.class,
                                                                 sample.cohort = values$sample.cohort,
                                                                 localT2_mode = input$I_perc_localT2,
                                                                 nWorkers = nWorkers)
        } else if (input$I_TH_METHOD=="singleTH") {
          values$O_txtout_I_GMCS_LIST <- gMCS.info$selected$fullname
          values$O_txtout_I_TH_METHOD <- paste0("singleTH (", input$I_perc_singleTH, ")")
          ResultsEssentiality <- CalculateEssentialGenes_singleTH(gene.exp = values$gene.exp,
                                                                  gMCS.info = gMCS.info$selected,
                                                                  sample.class = values$sample.class,
                                                                  singleTH = input$I_perc_singleTH,
                                                                  nWorkers = nWorkers)
        }
        # reorder results to have first genes more essential in all samples
        idx <- order(apply(ResultsEssentiality$ratio.essential.gene[,levels(values$sample.class)],1,mean), decreasing = T)  
        ResultsEssentiality$num.essential.gene <- ResultsEssentiality$num.essential.gene[idx,]
        ResultsEssentiality$ratio.essential.gene <- ResultsEssentiality$ratio.essential.gene[idx,]
        
        # eliminate rownames in tables, not informative
        rownames(ResultsEssentiality$num.essential.gene) <- 1:nrow(ResultsEssentiality$num.essential.gene)
        rownames(ResultsEssentiality$ratio.essential.gene) <- 1:nrow(ResultsEssentiality$ratio.essential.gene)
        
        # save parameters of calculation
        ResultsEssentiality$options <- list(
          I_GMCS_LIST = values$I_GMCS_LIST,
          I_GMCS_LIST_selected_tasks = metTasks[[values$I_GMCS_LIST]]$DESCRIPTION,
          I_GMCS_LIST_max_length = max(gMCS.info$selected$gMCSs.ENSEMBL.length),
          I_TH_METHOD = input$I_TH_METHOD,
          I_perc_gmcsTH = input$I_perc_gmcsTH,
          I_perc_localT2 = input$I_perc_localT2,
          I_perc_singleTH = input$I_perc_singleTH
        )
        
        # update and store results
        values$ResultsEssentiality <- ResultsEssentiality
        
        flags$flag_show_table_essentiality <- T
        flags$flag_show_table_gmcs_essentiality <- T
        flags$flag_show_table_gmcs_correlation <- T
        hideElement("mp3_loading")
      }
    })
    
    
    # decide which upload data to use ####
    showElement("mp3_ResultsEssentiality_input_RDATA")
    hideElement("mp3_ResultsEssentiality_input_example")
    
    observeEvent(input$I_ResultsEssentiality_show_input_examples,{
      hideElement("mp3_ResultsEssentiality_input_RDATA")
      showElement("mp3_ResultsEssentiality_input_example")
    })
    observeEvent(input$I_ResultsEssentiality_show_input_RDATA,{
      showElement("mp3_ResultsEssentiality_input_RDATA")
      hideElement("mp3_ResultsEssentiality_input_example")
    })
    
    
    # load precomputed results function ####
    load_precomputed_results_from_filename <- function(fname) {
      hide_all_elements(avoid = c("mp2_1", "mp2_2", "mp2_3", "mp2_4"))
      showElement("mp3_loading")
      flags$flag_show_summary_data <- F
      flags$flag_show_table_essentiality <- F
      flags$flag_show_heatmap_gmcs <- F
      flags$flag_show_dotplot_gmcs <- F
      
      load(fname)
      
      values$ResultsEssentiality <- ResultsEssentiality
      values$gene.exp  <- gene.exp
      values$sample.class <- sample.class
      values$sample.cohort <- sample.cohort
      
      idx <- order(apply(ResultsEssentiality$ratio.essential.gene[,levels(values$sample.class)],1,mean), decreasing = T)
      
      ResultsEssentiality$num.essential.gene <- ResultsEssentiality$num.essential.gene[idx,]
      ResultsEssentiality$ratio.essential.gene <- ResultsEssentiality$ratio.essential.gene[idx,]
      
      rownames(ResultsEssentiality$num.essential.gene) <- 1:nrow(ResultsEssentiality$num.essential.gene)
      rownames(ResultsEssentiality$ratio.essential.gene) <- 1:nrow(ResultsEssentiality$ratio.essential.gene)
      
      values$ResultsEssentiality <- ResultsEssentiality
      
      flags$flag_show_summary_data <- T
      flags$flag_show_table_essentiality <- T
      hideElement("mp3_loading")
      
      if (grepl('EssentialTasks',ResultsEssentiality$options$I_GMCS_LIST)) {
        updateRadioButtons(session, inputId = "I_GMCS_1_EssentialTask_or_Only", selected = 'EssentialTasks')
        updateRadioButtons(session, inputId = "I_GMCS_2_EssentialTask_or_Only", selected = 'EssentialTasks')
      } else if (grepl('Only',ResultsEssentiality$options$I_GMCS_LIST)) {
        updateRadioButtons(session, inputId = "I_GMCS_1_EssentialTask_or_Only", selected = 'Only')
        updateRadioButtons(session, inputId = "I_GMCS_2_EssentialTask_or_Only", selected = 'Only')
      } else {
        updateRadioButtons(session, inputId = "I_GMCS_1_EssentialTask_or_Only", selected = 'Custom')
        updateRadioButtons(session, inputId = "I_GMCS_2_EssentialTask_or_Only", selected = 'Custom')
        # change the sliders and the tree in the input according to the load tasks
      }
      if (grepl('Culture',ResultsEssentiality$options$I_GMCS_LIST)) {
        updateRadioButtons(session, inputId = "I_GMCS_1_CultureMedia_or_Full", selected = 'CultureMedium')
        updateRadioButtons(session, inputId = "I_GMCS_2_CultureMedia_or_Full", selected = 'CultureMedium')
      } else {
        updateRadioButtons(session, inputId = "I_GMCS_1_CultureMedia_or_Full", selected = 'FullMedium')
        updateRadioButtons(session, inputId = "I_GMCS_2_CultureMedia_or_Full", selected = 'FullMedium')
      }
      
      # updateRadioButtons(session, inputId = "I_GMCS_LIST", selected = ResultsEssentiality$options$I_GMCS_LIST)
      updateRadioButtons(session, inputId = "I_TH_METHOD", selected = ResultsEssentiality$options$I_TH_METHOD)
      updateSliderInput(session, inputId = "I_perc_gmcsTH", value = ResultsEssentiality$options$I_perc_gmcsTH)
      updateRadioButtons(session, inputId = "I_perc_localT2", selected = ResultsEssentiality$options$I_perc_localT2)
      updateNumericInput(session, inputId = "I_perc_singleTH", value = ResultsEssentiality$options$I_perc_singleTH)
      
      values$O_txtout_I_GMCS_LIST <- gMCS.info$selected$fullname
      if (ResultsEssentiality$options$I_TH_METHOD=="gmcsTH") {
        values$O_txtout_I_TH_METHOD <- paste0("gmcsTH (", input$I_perc_gmcsTH*100, "%)")
      } else if (ResultsEssentiality$options$I_TH_METHOD=="localT2") {
        values$O_txtout_I_TH_METHOD <- paste0("localT2 (", input$I_perc_localT2, ")")
      } else if (ResultsEssentiality$options$I_TH_METHOD=="singleTH") {
        values$O_txtout_I_TH_METHOD <- paste0("singleTH (", input$I_perc_singleTH, ")")
      }
      
      
      flags$flag_show_table_gmcs_essentiality <- T
      flags$flag_show_table_gmcs_correlation <- T
    }
    
    
    # load precomputed results options ####
    observeEvent(input$ResultsEssentiality_input,{
      if (is.null(input$ResultsEssentiality_input) == FALSE){
        print(input$ResultsEssentiality_input)
        load_precomputed_results_from_filename(input$ResultsEssentiality_input$datapath)
      }
    })
    
    # load precomputed results from example datasets [gMCST5]
    observeEvent(input$I_ResultsEssentiality_input_example_gmcsTH_TPM,{
      load_precomputed_results_from_filename("./ExampleData/ExampleResults_HumanGEM_BcellMM_gmcsTH5_TPM.Rdata")
      values$ResultsEssentiality$options <- list(I_GMCS_LIST = "EssentialTasks_CultureMedium",
                                                 I_GMCS_LIST_selected_tasks = metTasks$EssentialTasks_CultureMedium$DESCRIPTION,
                                                 I_GMCS_LIST_max_length = max(gMCS.info.raw$EssentialTasks_CultureMedium$gMCSs.ENSEMBL.length),
                                                 I_TH_METHOD = "gmcsTH",
                                                 I_perc_gmcsTH = 0.05,
                                                 I_perc_localT2 = "all_genes_gMCSs",
                                                 I_perc_singleTH = 1)
      # updateCheckboxInput(session, inputId = "I_isTarget_GEA_class_7", value = T)
      updateCheckboxGroupInput(session, inputId = "I_isTarget_GEA_class", choices = levels(values$sample.class), selected = c("MM"))
      flags$flag_show_table_gmcs_essentiality <- T
      flags$flag_show_table_gmcs_correlation <- T
      show_button_calculate_essentiality()
    })
    # load precomputed results from example datasets [localT2]
    observeEvent(input$I_ResultsEssentiality_input_example_localT2_TPM,{
      load_precomputed_results_from_filename("./ExampleData/ExampleResults_HumanGEM_BcellMM_localT2_TPM.Rdata")
      values$ResultsEssentiality$options <- list(I_GMCS_LIST = "EssentialTasks_CultureMedium",
                                                 I_GMCS_LIST_selected_tasks = metTasks$EssentialTasks_CultureMedium$DESCRIPTION,
                                                 I_GMCS_LIST_max_length = max(gMCS.info.raw$EssentialTasks_CultureMedium$gMCSs.ENSEMBL.length),
                                                 I_TH_METHOD = "localT2",
                                                 I_perc_gmcsTH = 0.05,
                                                 I_perc_localT2 = "all_genes_gMCSs",
                                                 I_perc_singleTH = 1)
      # updateCheckboxInput(session, inputId = "I_isTarget_GEA_class_7", value = T)
      updateCheckboxGroupInput(session, inputId = "I_isTarget_GEA_class", choices = levels(values$sample.class), selected = c("MM"))
      flags$flag_show_table_gmcs_essentiality <- T
      flags$flag_show_table_gmcs_correlation <- T
      show_button_calculate_essentiality()
    })
    # load precomputed results from all DepMap
    observeEvent(input$I_ResultsEssentiality_input_all_DepMap_gmcsTH_TPM,{
      load_precomputed_results_from_filename("./ExampleData/CompleteResults_HumanGEM_DepMap_gmcsTH5_TPM.Rdata")
      values$ResultsEssentiality$options <- list(I_GMCS_LIST = "EssentialTasks_CultureMedium",
                                                 I_GMCS_LIST_selected_tasks = metTasks$EssentialTasks_CultureMedium$DESCRIPTION,
                                                 I_GMCS_LIST_max_length = max(gMCS.info.raw$EssentialTasks_CultureMedium$gMCSs.ENSEMBL.length),
                                                 I_TH_METHOD = "gmcsTH",
                                                 I_perc_gmcsTH = 0.05,
                                                 I_perc_localT2 = "all_genes_gMCSs",
                                                 I_perc_singleTH = 1)
      flags$flag_show_table_gmcs_essentiality <- T
      flags$flag_show_table_gmcs_correlation <- T
      show_button_calculate_essentiality()
    })
    
    
    
    # load precomputed results from example datasets [gMCST5]
    observeEvent(input$I_ResultsEssentiality_input_example_gmcsTH_log2TPM,{
      load_precomputed_results_from_filename("./ExampleData/ExampleResults_HumanGEM_BcellMM_gmcsTH5_log2TPM.Rdata")
      values$ResultsEssentiality$options <- list(I_GMCS_LIST = "EssentialTasks_CultureMedium",
                                                 I_GMCS_LIST_selected_tasks = metTasks$EssentialTasks_CultureMedium$DESCRIPTION,
                                                 I_GMCS_LIST_max_length = max(gMCS.info.raw$EssentialTasks_CultureMedium$gMCSs.ENSEMBL.length),
                                                 I_TH_METHOD = "gmcsTH",
                                                 I_perc_gmcsTH = 0.05,
                                                 I_perc_localT2 = "all_genes_gMCSs",
                                                 I_perc_singleTH = 1)
      # updateCheckboxInput(session, inputId = "I_isTarget_GEA_class_7", value = T)
      updateCheckboxGroupInput(session, inputId = "I_isTarget_GEA_class", choices = levels(values$sample.class), selected = c("MM"))
      flags$flag_show_table_gmcs_essentiality <- T
      flags$flag_show_table_gmcs_correlation <- T
      show_button_calculate_essentiality()
    })
    # load precomputed results from example datasets [localT2]
    observeEvent(input$I_ResultsEssentiality_input_example_localT2_log2TPM,{
      load_precomputed_results_from_filename("./ExampleData/ExampleResults_HumanGEM_BcellMM_localT2_log2TPM.Rdata")
      values$ResultsEssentiality$options <- list(I_GMCS_LIST = "EssentialTasks_CultureMedium",
                                                 I_GMCS_LIST_selected_tasks = metTasks$EssentialTasks_CultureMedium$DESCRIPTION,
                                                 I_GMCS_LIST_max_length = max(gMCS.info.raw$EssentialTasks_CultureMedium$gMCSs.ENSEMBL.length),
                                                 I_TH_METHOD = "localT2",
                                                 I_perc_gmcsTH = 0.05,
                                                 I_perc_localT2 = "all_genes_gMCSs",
                                                 I_perc_singleTH = 1)
      # updateCheckboxInput(session, inputId = "I_isTarget_GEA_class_7", value = T)
      updateCheckboxGroupInput(session, inputId = "I_isTarget_GEA_class", choices = levels(values$sample.class), selected = c("MM"))
      flags$flag_show_table_gmcs_essentiality <- T
      flags$flag_show_table_gmcs_correlation <- T
      show_button_calculate_essentiality()
    })
    
    
    # calculate which table to show
    ResultsEssentialityTable <- reactive({
      if (flags$flag_show_table_essentiality && !is.null(values$ResultsEssentiality)){
        # browser()
        if (input$I_RESULT_TABLE_MODE=="percentage"){
          zzz <- values$ResultsEssentiality$ratio.essential.gene
          ResultsEssentialityTable <- datatable(zzz, filter = "top", rownames=F, selection = "single",
                                                options(pageLength = 10, lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "All"))))  %>%  formatPercentage(columns = levels(values$sample.class),digits=2)
        } else {
          zzz <- values$ResultsEssentiality$num.essential.gene
          ResultsEssentialityTable <- datatable(zzz, filter = "top", rownames=F, selection = "single",
                                                options(pageLength = 10, lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "All"))))
        }
        # rownames(ResultsEssentialityTable) <- NULL
      } else {
        hideElement("mp3")
        ResultsEssentialityTable <- NULL
      }
      return(ResultsEssentialityTable)
    })
    
    # calculate which table to show
    ResultsEssentialityTableToTxt <- reactive({
      if (flags$flag_show_table_essentiality && !is.null(values$ResultsEssentiality)){
        if (input$I_RESULT_TABLE_MODE=="percentage"){
          ResultsEssentialityTableToTxt <- values$ResultsEssentiality$ratio.essential.gene
        } else {
          ResultsEssentialityTableToTxt <- values$ResultsEssentiality$num.essential.gene
        }
      } else {
        ResultsEssentialityTableToTxt <- NULL
      }
    })
    
    
    # browser()
    # Table of essential genes
    observeEvent(flags$flag_show_table_essentiality, {
      print(paste("flags$flag_show_table_essentiality = ", flags$flag_show_table_essentiality))
      if (!flags$flag_show_table_essentiality){
        hideElement("mp3")
        hideElement("mp3_loading")
        output$O_table_essential_genes <- NULL
      } else {
        showElement("mp3")
        output$O_table_essential_genes <- DT::renderDataTable({ResultsEssentialityTable()}) 
      }
    })
    
    # Save the results ####
    output$O_results_essential_genes_simple_download <- downloadHandler(filename = paste0("gMCStool_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".csv"), 
                                                                        content = function(fname){fwrite(ResultsEssentialityTableToTxt(), file = fname, quote = F, row.names = F, sep = ";")})
    output$O_results_essential_genes_all_download <- downloadHandler(filename = paste0("gMCStool_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".xlsx"), 
                                                                     content = function(fname){SaveResultsInExcel(values$ResultsEssentiality, gMCS.info$selected, filename = fname)})
    output$O_results_all_download <- downloadHandler(filename = paste0("gMCStool_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".Rdata"), 
                                                     content = function(fname){ResultsEssentiality <- values$ResultsEssentiality;gene.exp <- values$gene.exp; sample.class <- values$sample.class; sample.cohort <- values$sample.cohort;
                                                     save(ResultsEssentiality, gene.exp, sample.class, sample.cohort, file = fname)})
    
    
    observeEvent(input$O_table_essential_genes_rows_selected, {
      print(paste0("Gene table row: ", input$O_table_essential_genes_rows_selected))
    })
  }
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # 4: CASE BY CASE ####
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  {
    # Define auxiliary functions to hide table and content
    show_button_calculate_essentiality <- function(){
      if (!RealTimeTables_mp4){
        hideElement("mp4")
        flags$flag_show_table_gmcs_essentiality <- F
        flags$flag_show_heatmap_gmcs <- F
        showElement("mp4_button")
      } else {
        hideElement("mp4_button")
        showElement("mp4")
      }
    }
    
    # Code to show the table of correlations
    observeEvent(input$I_mp4_calculate, {
      flags$flag_show_table_gmcs_essentiality <- T
      flags$flag_show_table_gmcs_correlation <- T
      hideElement("mp4_button")
      showElement("mp4")
    })
    
    
    
    # Define the sliders and the table to show ####
    observeEvent(values$sample.class, {
      for(i in 1:num_max_classes) { hideElement(paste0("mp4_class_",i)) }
      
      for(i in 1:length(levels(values$sample.class))) { 
        showElement(paste0("mp4_class_",i))
      }
      if (length(levels(values$sample.class))>=1) output$O_name_GEA_class_1 <- renderText({levels(values$sample.class)[1]})
      if (length(levels(values$sample.class))>=2) output$O_name_GEA_class_2 <- renderText({levels(values$sample.class)[2]})
      if (length(levels(values$sample.class))>=3) output$O_name_GEA_class_3 <- renderText({levels(values$sample.class)[3]})
      if (length(levels(values$sample.class))>=4) output$O_name_GEA_class_4 <- renderText({levels(values$sample.class)[4]})
      if (length(levels(values$sample.class))>=5) output$O_name_GEA_class_5 <- renderText({levels(values$sample.class)[5]})
      if (length(levels(values$sample.class))>=6) output$O_name_GEA_class_6 <- renderText({levels(values$sample.class)[6]})
      if (length(levels(values$sample.class))>=7) output$O_name_GEA_class_7 <- renderText({levels(values$sample.class)[7]})
      if (length(levels(values$sample.class))>=8) output$O_name_GEA_class_8 <- renderText({levels(values$sample.class)[8]})
      if (length(levels(values$sample.class))>=9) output$O_name_GEA_class_9 <- renderText({levels(values$sample.class)[9]})
      if (length(levels(values$sample.class))>=10) output$O_name_GEA_class_10 <- renderText({levels(values$sample.class)[10]})
      if (length(levels(values$sample.class))>=11) output$O_name_GEA_class_11 <- renderText({levels(values$sample.class)[11]})
      if (length(levels(values$sample.class))>=12) output$O_name_GEA_class_12 <- renderText({levels(values$sample.class)[12]})
      if (length(levels(values$sample.class))>=13) output$O_name_GEA_class_13 <- renderText({levels(values$sample.class)[13]})
      if (length(levels(values$sample.class))>=14) output$O_name_GEA_class_14 <- renderText({levels(values$sample.class)[14]})
      if (length(levels(values$sample.class))>=15) output$O_name_GEA_class_15 <- renderText({levels(values$sample.class)[15]})
      if (length(levels(values$sample.class))>=16) output$O_name_GEA_class_16 <- renderText({levels(values$sample.class)[16]})
      if (length(levels(values$sample.class))>=17) output$O_name_GEA_class_17 <- renderText({levels(values$sample.class)[17]})
      if (length(levels(values$sample.class))>=18) output$O_name_GEA_class_18 <- renderText({levels(values$sample.class)[18]})
      if (length(levels(values$sample.class))>=19) output$O_name_GEA_class_19 <- renderText({levels(values$sample.class)[19]})
      if (length(levels(values$sample.class))>=20) output$O_name_GEA_class_20 <- renderText({levels(values$sample.class)[20]})
      if (length(levels(values$sample.class))>=21) output$O_name_GEA_class_21 <- renderText({levels(values$sample.class)[21]})
      if (length(levels(values$sample.class))>=22) output$O_name_GEA_class_22 <- renderText({levels(values$sample.class)[22]})
      if (length(levels(values$sample.class))>=23) output$O_name_GEA_class_23 <- renderText({levels(values$sample.class)[23]})
      if (length(levels(values$sample.class))>=24) output$O_name_GEA_class_24 <- renderText({levels(values$sample.class)[24]})
      if (length(levels(values$sample.class))>=25) output$O_name_GEA_class_25 <- renderText({levels(values$sample.class)[25]})
      if (length(levels(values$sample.class))>=26) output$O_name_GEA_class_26 <- renderText({levels(values$sample.class)[26]})
      if (length(levels(values$sample.class))>=27) output$O_name_GEA_class_27 <- renderText({levels(values$sample.class)[27]})
      if (length(levels(values$sample.class))>=28) output$O_name_GEA_class_28 <- renderText({levels(values$sample.class)[28]})
      if (length(levels(values$sample.class))>=29) output$O_name_GEA_class_29 <- renderText({levels(values$sample.class)[29]})
      if (length(levels(values$sample.class))>=30) output$O_name_GEA_class_30 <- renderText({levels(values$sample.class)[30]})
      
      # update the possible target class
      updateCheckboxGroupInput(session, inputId = "I_isTarget_GEA_class", 
                               choices = levels(values$sample.class),
                               selected = NULL)
    })
    
    
    observeEvent(input$I_RESULT_TABLE_GMCS_MODE, {
      updateRadioButtons(session, inputId = "I_RESULT_TABLE_GMCS_MODE_bis", selected = input$I_RESULT_TABLE_GMCS_MODE)
      flags$flag_show_heatmap_gmcs <- F
      flags$flag_show_dotplot_gmcs <- F
      if (input$I_RESULT_TABLE_GMCS_MODE=="percentage"){
        for (ii in 1:length(levels(values$sample.class))) { 
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii), min = 0, max = 100, value = c(0,100), step = 1)
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii, "_bis"), min = 0, max = 100, value = c(0,100), step = 1)
        }
        updateSliderInput(session, "I_perc_GEA_class_non_target_filter", min = 0, max = 100, value = c(0,0), step = 1)
        updateSliderInput(session, "I_perc_GEA_class_target_filter", min = 0, max = 100, value = c(1,100), step = 1)
      } else {
        for (ii in 1:length(levels(values$sample.class))) {
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii), min = 0,
                            max = sum(values$sample.class == levels(values$sample.class)[ii]),
                            value = c(0,sum(values$sample.class == levels(values$sample.class)[ii])), step = 1)
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii,"_bis"), min = 0,
                            max = sum(values$sample.class == levels(values$sample.class)[ii]),
                            value = c(0,sum(values$sample.class == levels(values$sample.class)[ii])), step = 1)
        }
        # browser()
        updateSliderInput(session, "I_perc_GEA_class_non_target_filter", min = 0, max = min(table(values$sample.class)[!values$sample.class.target[1:length(levels(values$sample.class))]]), value = c(0,0), step = 1)
        updateSliderInput(session, "I_perc_GEA_class_target_filter", min = 0, max = min(table(values$sample.class)[values$sample.class.target]), value = c(1,100), step = 1)
      }
    })
    
    # Update the table ####
    observeEvent(c(input$I_perc_GEA_class_1, input$I_perc_GEA_class_2,
                   input$I_perc_GEA_class_3, input$I_perc_GEA_class_4,
                   input$I_perc_GEA_class_5, input$I_perc_GEA_class_6,
                   input$I_perc_GEA_class_7, input$I_perc_GEA_class_8,
                   input$I_perc_GEA_class_9, input$I_perc_GEA_class_10,
                   input$I_perc_GEA_class_11, input$I_perc_GEA_class_12,
                   input$I_perc_GEA_class_13, input$I_perc_GEA_class_14,
                   input$I_perc_GEA_class_15, input$I_perc_GEA_class_16,
                   input$I_perc_GEA_class_17, input$I_perc_GEA_class_18,
                   input$I_perc_GEA_class_19, input$I_perc_GEA_class_20,
                   input$I_perc_GEA_class_21, input$I_perc_GEA_class_22,
                   input$I_perc_GEA_class_23, input$I_perc_GEA_class_24,
                   input$I_perc_GEA_class_25, input$I_perc_GEA_class_26,
                   input$I_perc_GEA_class_27, input$I_perc_GEA_class_28,
                   input$I_perc_GEA_class_29, input$I_perc_GEA_class_30),{
                     
                     # store results if they exist
                     if (length(levels(values$sample.class))>=1) values$filtersGEA[[1]] <- input$I_perc_GEA_class_1
                     if (length(levels(values$sample.class))>=2) values$filtersGEA[[2]] <- input$I_perc_GEA_class_2
                     if (length(levels(values$sample.class))>=3) values$filtersGEA[[3]] <- input$I_perc_GEA_class_3
                     if (length(levels(values$sample.class))>=4) values$filtersGEA[[4]] <- input$I_perc_GEA_class_4
                     if (length(levels(values$sample.class))>=5) values$filtersGEA[[5]] <- input$I_perc_GEA_class_5
                     if (length(levels(values$sample.class))>=6) values$filtersGEA[[6]] <- input$I_perc_GEA_class_6
                     if (length(levels(values$sample.class))>=7) values$filtersGEA[[7]] <- input$I_perc_GEA_class_7
                     if (length(levels(values$sample.class))>=8) values$filtersGEA[[8]] <- input$I_perc_GEA_class_8
                     if (length(levels(values$sample.class))>=9) values$filtersGEA[[9]] <- input$I_perc_GEA_class_9
                     if (length(levels(values$sample.class))>=10) values$filtersGEA[[10]] <- input$I_perc_GEA_class_10
                     if (length(levels(values$sample.class))>=11) values$filtersGEA[[11]] <- input$I_perc_GEA_class_11
                     if (length(levels(values$sample.class))>=12) values$filtersGEA[[12]] <- input$I_perc_GEA_class_12
                     if (length(levels(values$sample.class))>=13) values$filtersGEA[[13]] <- input$I_perc_GEA_class_13
                     if (length(levels(values$sample.class))>=14) values$filtersGEA[[14]] <- input$I_perc_GEA_class_14
                     if (length(levels(values$sample.class))>=15) values$filtersGEA[[15]] <- input$I_perc_GEA_class_15
                     if (length(levels(values$sample.class))>=16) values$filtersGEA[[16]] <- input$I_perc_GEA_class_16
                     if (length(levels(values$sample.class))>=17) values$filtersGEA[[17]] <- input$I_perc_GEA_class_17
                     if (length(levels(values$sample.class))>=18) values$filtersGEA[[18]] <- input$I_perc_GEA_class_18
                     if (length(levels(values$sample.class))>=19) values$filtersGEA[[19]] <- input$I_perc_GEA_class_19
                     if (length(levels(values$sample.class))>=20) values$filtersGEA[[20]] <- input$I_perc_GEA_class_20
                     if (length(levels(values$sample.class))>=21) values$filtersGEA[[21]] <- input$I_perc_GEA_class_21
                     if (length(levels(values$sample.class))>=22) values$filtersGEA[[22]] <- input$I_perc_GEA_class_22
                     if (length(levels(values$sample.class))>=23) values$filtersGEA[[23]] <- input$I_perc_GEA_class_23
                     if (length(levels(values$sample.class))>=24) values$filtersGEA[[24]] <- input$I_perc_GEA_class_24
                     if (length(levels(values$sample.class))>=25) values$filtersGEA[[25]] <- input$I_perc_GEA_class_25
                     if (length(levels(values$sample.class))>=26) values$filtersGEA[[26]] <- input$I_perc_GEA_class_26
                     if (length(levels(values$sample.class))>=27) values$filtersGEA[[27]] <- input$I_perc_GEA_class_27
                     if (length(levels(values$sample.class))>=28) values$filtersGEA[[28]] <- input$I_perc_GEA_class_28
                     if (length(levels(values$sample.class))>=29) values$filtersGEA[[29]] <- input$I_perc_GEA_class_29
                     if (length(levels(values$sample.class))>=30) values$filtersGEA[[30]] <- input$I_perc_GEA_class_30
                     
                     
                     # as the tables are going to change, remove selection of table
                     values$O_table_essential_gmcs_bis_rows_selected_previous_selection <- 0
                     values$O_table_essential_gmcs_rows_selected_previous_selection <- 0
                     
                     # update sliders of the other page
                     for (ii in 1:length(levels(values$sample.class))) {
                       updateSliderInput(session, paste0("I_perc_GEA_class_",ii, "_bis"), value = values$filtersGEA[[ii]])
                     }
                     
                     # # if not real time, hide tables and show calculation button
                     show_button_calculate_essentiality()
                   })
    
    
    # Select classes that are essential
    observeEvent(c(input$I_isTarget_GEA_class),{
      # store the new values
      values$sample.class.target <- levels(values$sample.class) %in% input$I_isTarget_GEA_class
      
      # update the general sliders
      if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
        # browser()
        updateSliderInput(session, "I_perc_GEA_class_non_target_filter", min = 0, max = min(table(values$sample.class)[!values$sample.class.target[1:length(levels(values$sample.class))]]), value = c(0,0), step = 1)
        updateSliderInput(session, "I_perc_GEA_class_target_filter", min = 0, max = min(table(values$sample.class)[values$sample.class.target]), value = c(1,100), step = 1)
      }
      
    })
    
    # udate sliders according to filters
    observeEvent(input$I_class_1_30_filter_targeted_class, {
      for (ii in 1:length(levels(values$sample.class))) {
        if (values$sample.class.target[ii]==0) {
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii), value = input$I_perc_GEA_class_non_target_filter)
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii, "_bis"), value = input$I_perc_GEA_class_non_target_filter)
        } else {
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii), value = input$I_perc_GEA_class_target_filter)
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii, "_bis"), value = input$I_perc_GEA_class_target_filter)
        }
      }
    })
    
    # reset the sliders
    observeEvent(input$I_class_1_30_unfilter_targeted_class, {
      if (input$I_RESULT_TABLE_GMCS_MODE=="percentage"){
        for (ii in 1:length(levels(values$sample.class))) { 
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii), min = 0, max = 100, value = c(0,100), step = 1)
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii, "_bis"), min = 0, max = 100, value = c(0,100), step = 1)
        }
      } else {
        for (ii in 1:length(levels(values$sample.class))) {
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii), min = 0,
                            max = sum(values$sample.class == levels(values$sample.class)[ii]),
                            value = c(0,sum(values$sample.class == levels(values$sample.class)[ii])), step = 1)
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii,"_bis"), min = 0,
                            max = sum(values$sample.class == levels(values$sample.class)[ii]),
                            value = c(0,sum(values$sample.class == levels(values$sample.class)[ii])), step = 1)
        }
      }
    })
    
    
    # calculate which table to show, taking into account all gMCS at a time for the filters ####
    ResultsGMCSsTable <- reactive({
      if (flags$flag_show_table_gmcs_essentiality){
        list.gMCS.essential.aux <- values$ResultsEssentiality$list.gMCS.essential
        # preprocess to obtain ratio or number
        if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
          list.gene.essential.aux <- values$ResultsEssentiality$num.essential.gene
        } else {
          list.gene.essential.aux <- values$ResultsEssentiality$ratio.essential.gene
        }
        
        if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
          for (i in 1:length(levels(values$sample.class)))
          {
            list.gMCS.essential.aux[,levels(values$sample.class)[i]] = list.gMCS.essential.aux[,levels(values$sample.class)[i]] * sum(values$sample.class == levels(values$sample.class)[i])
          }
        }
        for (i in 1:length(levels(values$sample.class)))
        {
          if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
            idx <- (list.gene.essential.aux[,levels(values$sample.class)[i]] >= values$filtersGEA[[i]][1]) &
              (list.gene.essential.aux[,levels(values$sample.class)[i]] <= values$filtersGEA[[i]][2])
          } else {
            idx <- (list.gene.essential.aux[,levels(values$sample.class)[i]] >= values$filtersGEA[[i]][1]/100) &
              (list.gene.essential.aux[,levels(values$sample.class)[i]] <= values$filtersGEA[[i]][2]/100)
          }
          list.gene.essential.aux <- list.gene.essential.aux[idx,]
        }
        list.gMCS.essential.aux <- list.gMCS.essential.aux %>% filter(ENSEMBL %in% list.gene.essential.aux$ENSEMBL)
        rownames(list.gMCS.essential.aux) <- 1:nrow(list.gMCS.essential.aux)
        print(head(rownames(list.gMCS.essential.aux)))
        values$list.gMCS.essential.filtered <- list.gMCS.essential.aux
        if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
          ResultsGMCSsTable <- datatable(values$list.gMCS.essential.filtered, filter = "top", selection = "single", rownames = F,
                                         options(pageLength = 10, lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "All"))))
        } else {
          ResultsGMCSsTable <- datatable(values$list.gMCS.essential.filtered, filter = "top", selection = "single", rownames = F,
                                         options(pageLength = 10, lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "All"))))   %>%  formatPercentage(columns = levels(values$sample.class),digits=2)
        }
      } else {
        ResultsGMCSsTable <- NULL
      }
    })
    
    # Show gMCS essentials
    observeEvent(flags$flag_show_table_gmcs_essentiality, {
      if (!flags$flag_show_table_gmcs_essentiality){
        print(paste("flags$flag_show_table_essentiality =", flags$flag_show_table_essentiality))
        print(paste("flags$flag_show_heatmap_gmcs =", flags$flag_show_heatmap_gmcs))
        hideElement("mp4"); hideElement("mp4_heatmap_1"); hideElement("mp4_heatmap_2"); hideElement("mp4_heatmap_3") 
        hideElement("mp5"); hideElement("mp5_boxplot")
        # showElement("mp5")
        output$O_table_essential_gmcs <- NULL
      } else {
        showElement("mp4")
        # showElement("mp5")
        output$O_table_essential_gmcs <- DT::renderDataTable({ResultsGMCSsTable()}
        )
      }
    })
    
    observeEvent(input$O_table_essential_gmcs_rows_selected, {
      # browser()
      if (values$O_table_essential_gmcs_rows_selected_previous_selection!=input$O_table_essential_gmcs_rows_selected) {
        
        flags$flag_show_heatmap_gmcs <- F
        flags$flag_show_dotplot_gmcs <- F
        print(paste0("gMCSs table row: ", input$O_table_essential_gmcs_rows_selected))
        print(values$list.gMCS.essential.filtered[input$O_table_essential_gmcs_rows_selected,])
        if (!is.null(values$list.gMCS.essential.filtered)){ flags$flag_show_heatmap_gmcs <- T }
        values$O_table_essential_gmcs_rows_selected_previous_selection=input$O_table_essential_gmcs_rows_selected
        
        flags$flag_show_heatmap_gmcs <- T
        flags$flag_show_dotplot_gmcs <- F
      } else {
        print(paste0("gMCSs table row: ", input$O_table_essential_gmcs_rows_selected, ", same as before"))
        
      }
    })
    
    # show the selected gMCS for certain gene ####
    observeEvent(c(input$I_include_boxplots_heatmap,
                   input$I_show_SYMBOL_heatmap,
                   input$I_colors_heatmap,
                   input$I_colors_colors_annotation,
                   input$I_isTarget_GEA_class_1, input$I_isTarget_GEA_class_2,
                   input$I_isTarget_GEA_class_3, input$I_isTarget_GEA_class_4,
                   input$I_isTarget_GEA_class_5, input$I_isTarget_GEA_class_6,
                   input$I_isTarget_GEA_class_7, input$I_isTarget_GEA_class_8,
                   input$I_isTarget_GEA_class_9, input$I_isTarget_GEA_class_10,
                   input$I_isTarget_GEA_class_11, input$I_isTarget_GEA_class_12,
                   input$I_isTarget_GEA_class_13, input$I_isTarget_GEA_class_14,
                   input$I_isTarget_GEA_class_15, 
                   flags$flag_show_heatmap_gmcs), {
                     
                     # browser()
                     if (flags$flag_show_heatmap_gmcs){
                       # showElement("mp4_heatmap")
                       showElement("mp4_loading"); hideElement("mp4_heatmap"); hideElement("mp4_heatmap_1"); hideElement("mp4_heatmap_2")
                       output$O_table_essential_gmcs_selected_tasks <- renderTable({ResultsGMCSsTable_task()}, align = "l", bordered = T)
                       output$O_table_essential_gmcs_selected_mets <- renderTable({ResultsGMCSsTable_mets()}, align = "l", bordered = T)
                       # values$plot.obj <- show_heatmap()
                       if (input$I_TH_METHOD %in% c("gmcsTH", "singleTH")){
                         print("start ShowHeatmap_gmcsTH")
                         values$plot.obj <- ShowHeatmap_gmcsTH(values$gene.exp,
                                                               values$ResultsEssentiality$gene.ratio,
                                                               gMCS.info$selected,
                                                               values$list.gMCS.essential.filtered[input$O_table_essential_gmcs_rows_selected,],
                                                               values$sample.class,
                                                               values$sample.class.target,
                                                               values$sample.cohort,
                                                               input$I_include_boxplots_heatmap,
                                                               input$I_show_SYMBOL_heatmap,
                                                               input$I_colors_heatmap,
                                                               input$I_colors_colors_annotation,
                                                               input$I_RESULT_TABLE_GMCS_MODE)
                         print("end ShowHeatmap_gmcsTH")
                       } else {
                         print("start ShowHeatmap_localT2")
                         values$plot.obj <- ShowHeatmap_localT2(values$gene.exp,
                                                                values$ResultsEssentiality$localT2,
                                                                gMCS.info$selected,
                                                                values$list.gMCS.essential.filtered[input$O_table_essential_gmcs_rows_selected,],
                                                                values$sample.class,
                                                                values$sample.class.target,
                                                                values$sample.cohort,
                                                                input$I_include_boxplots_heatmap,
                                                                input$I_show_SYMBOL_heatmap,
                                                                input$I_colors_heatmap,
                                                                input$I_colors_colors_annotation,
                                                                input$I_RESULT_TABLE_GMCS_MODE)
                         print("end ShowHeatmap_localT2")
                       }
                       if (input$I_include_boxplots_heatmap){
                         hideElement("mp4_loading")
                         showElement("mp4_heatmap")
                         showElement("mp4_heatmap_1")
                         showElement("mp4_heatmap_2")
                         showElement("mp4_heatmap_3")
                         output$O_heatmap_seleted_gmcs <- renderPlot({values$plot.obj$heatmap})
                         output$O_boxplot_seleted_gmcs <- renderPlot({values$plot.obj$boxplot})
                       } else {
                         hideElement("mp4_loading")
                         showElement("mp4_heatmap")
                         showElement("mp4_heatmap_1")
                         hideElement("mp4_heatmap_2")
                         showElement("mp4_heatmap_3")
                         output$O_heatmap_seleted_gmcs <- renderPlot({values$plot.obj$heatmap})
                         output$O_boxplot_seleted_gmcs <- renderPlot({ggplot() + theme_void()})
                       }
                     } else {
                       hideElement("mp4_loading")
                       hideElement("mp4_heatmap")
                       hideElement("mp4_heatmap_1")
                       hideElement("mp4_heatmap_2")
                       hideElement("mp4_heatmap_3")
                       values$plot.obj <- list(ggplot() + theme_void(), ggplot() + theme_void())
                       output$O_heatmap_seleted_gmcs <- renderPlot({ggplot() + theme_void()})
                       output$O_boxplot_seleted_gmcs <- renderPlot({ggplot() + theme_void()}) 
                     }
                   })
    
    
    
    
    # Show the metabolic task and the metabolites inhibited by this gMCS
    ResultsGMCSsTable_task <- reactive({
      if (flags$flag_show_heatmap_gmcs){
        
        if (input$O_table_essential_gmcs_rows_selected != values$O_table_essential_gmcs_rows_selected_previous_selection){
          index <- values$list.gMCS.essential.filtered[input$O_table_essential_gmcs_rows_selected,]
          index <- as.numeric(as.character(index$gMCS))
        } else {
          index <- values$list.gMCS.essential.filtered[values$O_table_essential_gmcs_rows_selected_previous_selection,]
          index <- as.numeric(as.character(index$gMCS))
        }
        
        print(paste("index =",index))
        print(paste("O_table_essential_gmcs_rows_selected =",input$O_table_essential_gmcs_rows_selected))
        table.gMCSs <- gMCS.info$selected$table.gMCSs
        
        idx <- table.gMCSs$task[table.gMCSs$idx==index]
        
        print(paste("idx =",idx))
        
        ResultsGMCSsTable_task <- data.frame("task number" = metTasks[[values$I_GMCS_LIST]][idx,1],
                                             "gMCS index in task" = table.gMCSs$gMCS[table.gMCSs$idx==index],
                                             "task name" = metTasks[[values$I_GMCS_LIST]][idx,2],
                                             "task group" = metTasks[[values$I_GMCS_LIST]][idx,3])
        print(str(ResultsGMCSsTable_task))
        ResultsGMCSsTable_task <- as.data.frame(ResultsGMCSsTable_task)
      } else {
        ResultsGMCSsTable_task <- NULL
      }
    })
    
    ResultsGMCSsTable_mets <- reactive({
      if (flags$flag_show_heatmap_gmcs){
        
        if (input$O_table_essential_gmcs_rows_selected != values$O_table_essential_gmcs_rows_selected_previous_selection){
          index <- values$list.gMCS.essential.filtered[input$O_table_essential_gmcs_rows_selected,]
          gene <- as.character(index$ENSEMBL)
          index <- as.numeric(as.character(index$gMCS))
        } else {
          index <- values$list.gMCS.essential.filtered[values$O_table_essential_gmcs_rows_selected_previous_selection,]
          gene <- as.character(index$ENSEMBL)
          index <- as.numeric(as.character(index$gMCS))
        }
        
        # print(str(index))
        
        if (grepl('Culture', values$I_GMCS_LIST)){
          ResultsGMCSsTable_mets <- metsBiomass$CultureMedium
        } else {
          ResultsGMCSsTable_mets <- metsBiomass$FullMedium
        }
        
        # browser()
        ResultsGMCSsTable_mets <- ResultsGMCSsTable_mets %>% dplyr::filter(gMCSText==gMCS.info$selected$gMCSs.ENSEMBL.txt[index])
        
        if (nrow(ResultsGMCSsTable_mets)>0){
          ResultsGMCSsTable_mets <- data.frame("mets" = ResultsGMCSsTable_mets$mets,
                                               "metNames" = ResultsGMCSsTable_mets$metNames)
          ResultsGMCSsTable_mets <- as.data.frame(ResultsGMCSsTable_mets)
        } else {
          ResultsGMCSsTable_mets <- NULL
        }
        
        print(str(ResultsGMCSsTable_mets))
        ResultsGMCSsTable_mets
      } else {
        ResultsGMCSsTable_mets <- NULL
      }
    })
    
    
    
    # Save the results ####
    output$O_heatmap_seleted_gmcs_download_1 <- downloadHandler(filename = paste0("gMCStool_heatmap_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".png"), 
                                                                content = function(fname){ggsave(ggarrange(plotlist = values$plot.obj, nrow = length(values$plot.obj), ncol = 1), 
                                                                                                 file = fname,
                                                                                                 width = 10, height = ifelse(length(values$plot.obj)==1, 6, 9),
                                                                                                 units = "in", dpi = 400, bg = "white")})
    output$O_heatmap_seleted_gmcs_download_2 <- downloadHandler(filename = paste0("gMCStool_heatmap_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".jpg"), 
                                                                content = function(fname){ggsave(ggarrange(plotlist = values$plot.obj, nrow = length(values$plot.obj), ncol = 1), 
                                                                                                 file = fname,
                                                                                                 width = 10, height = ifelse(length(values$plot.obj)==1, 6, 9),
                                                                                                 units = "in", dpi = 400, bg = "white")})
    output$O_heatmap_seleted_gmcs_download_3 <- downloadHandler(filename = paste0("gMCStool_heatmap_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".RDS"), 
                                                                content = function(fname){saveRDS(values$plot.obj, file = fname)})
    
  }
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # 5: CASE BY CASE  Achilles####
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  {
    output$O_txtout_I_GMCS_LIST <- renderText({values$O_txtout_I_GMCS_LIST})
    output$O_txtout_I_GMCS_LIST_bis <- renderText({values$O_txtout_I_GMCS_LIST})
    output$O_txtout_I_TH_METHOD <- renderText({values$O_txtout_I_TH_METHOD})
    output$O_txtout_I_TH_METHOD_bis <- renderText({values$O_txtout_I_TH_METHOD})
    
    
    # Define auxiliary functions to hide table and content
    show_button_calculate_correlations <- function(){
      if (!RealTimeTables_mp5){
        hideElement("mp5")
        flags$flag_show_table_gmcs_correlation <- F
        flags$flag_show_dotplot_gmcs <- F
        showElement("mp5_button")
      } else {
        hideElement("mp5_button")
        showElement("mp5")
      }
    }
    
    # Code to show the table of correlations
    observeEvent(input$I_mp5_calculate, {
      flags$flag_show_table_gmcs_essentiality <- T
      flags$flag_show_table_gmcs_correlation <- T
      filter_DepMap_tables()
      hideElement("mp5_button")
      showElement("mp5")
    })
    
    
    # Define the sliders and the table to show ####
    observeEvent(values$sample.class, {
      for(i in 1:num_max_classes)  { hideElement(paste0("mp5_class_",i)) }
      # browser()
      for(i in 1:length(levels(values$sample.class))) { showElement(paste0("mp5_class_",i)) }
      
      if (length(levels(values$sample.class))>=1) output$O_name_GEA_class_1_bis <- renderText({levels(values$sample.class)[1]})
      if (length(levels(values$sample.class))>=2) output$O_name_GEA_class_2_bis <- renderText({levels(values$sample.class)[2]})
      if (length(levels(values$sample.class))>=3) output$O_name_GEA_class_3_bis <- renderText({levels(values$sample.class)[3]})
      if (length(levels(values$sample.class))>=4) output$O_name_GEA_class_4_bis <- renderText({levels(values$sample.class)[4]})
      if (length(levels(values$sample.class))>=5) output$O_name_GEA_class_5_bis <- renderText({levels(values$sample.class)[5]})
      if (length(levels(values$sample.class))>=6) output$O_name_GEA_class_6_bis <- renderText({levels(values$sample.class)[6]})
      if (length(levels(values$sample.class))>=7) output$O_name_GEA_class_7_bis <- renderText({levels(values$sample.class)[7]})
      if (length(levels(values$sample.class))>=8) output$O_name_GEA_class_8_bis <- renderText({levels(values$sample.class)[8]})
      if (length(levels(values$sample.class))>=9) output$O_name_GEA_class_9_bis <- renderText({levels(values$sample.class)[9]})
      if (length(levels(values$sample.class))>=10) output$O_name_GEA_class_10_bis <- renderText({levels(values$sample.class)[10]})
      if (length(levels(values$sample.class))>=11) output$O_name_GEA_class_11_bis <- renderText({levels(values$sample.class)[11]})
      if (length(levels(values$sample.class))>=12) output$O_name_GEA_class_12_bis <- renderText({levels(values$sample.class)[12]})
      if (length(levels(values$sample.class))>=13) output$O_name_GEA_class_13_bis <- renderText({levels(values$sample.class)[13]})
      if (length(levels(values$sample.class))>=14) output$O_name_GEA_class_14_bis <- renderText({levels(values$sample.class)[14]})
      if (length(levels(values$sample.class))>=15) output$O_name_GEA_class_15_bis <- renderText({levels(values$sample.class)[15]})
      if (length(levels(values$sample.class))>=16) output$O_name_GEA_class_16_bis <- renderText({levels(values$sample.class)[16]})
      if (length(levels(values$sample.class))>=17) output$O_name_GEA_class_17_bis <- renderText({levels(values$sample.class)[17]})
      if (length(levels(values$sample.class))>=18) output$O_name_GEA_class_18_bis <- renderText({levels(values$sample.class)[18]})
      if (length(levels(values$sample.class))>=19) output$O_name_GEA_class_19_bis <- renderText({levels(values$sample.class)[19]})
      if (length(levels(values$sample.class))>=20) output$O_name_GEA_class_20_bis <- renderText({levels(values$sample.class)[20]})
      if (length(levels(values$sample.class))>=21) output$O_name_GEA_class_21_bis <- renderText({levels(values$sample.class)[21]})
      if (length(levels(values$sample.class))>=22) output$O_name_GEA_class_22_bis <- renderText({levels(values$sample.class)[22]})
      if (length(levels(values$sample.class))>=23) output$O_name_GEA_class_23_bis <- renderText({levels(values$sample.class)[23]})
      if (length(levels(values$sample.class))>=24) output$O_name_GEA_class_24_bis <- renderText({levels(values$sample.class)[24]})
      if (length(levels(values$sample.class))>=25) output$O_name_GEA_class_25_bis <- renderText({levels(values$sample.class)[25]})
      if (length(levels(values$sample.class))>=26) output$O_name_GEA_class_26_bis <- renderText({levels(values$sample.class)[26]})
      if (length(levels(values$sample.class))>=27) output$O_name_GEA_class_27_bis <- renderText({levels(values$sample.class)[27]})
      if (length(levels(values$sample.class))>=28) output$O_name_GEA_class_28_bis <- renderText({levels(values$sample.class)[28]})
      if (length(levels(values$sample.class))>=29) output$O_name_GEA_class_29_bis <- renderText({levels(values$sample.class)[29]})
      if (length(levels(values$sample.class))>=30) output$O_name_GEA_class_30_bis <- renderText({levels(values$sample.class)[30]})
    })
    
    
    observeEvent(input$I_RESULT_TABLE_GMCS_MODE_bis, {
      updateRadioButtons(session = session, "I_RESULT_TABLE_GMCS_MODE", selected = input$I_RESULT_TABLE_GMCS_MODE_bis)
      flags$flag_show_heatmap_gmcs <- F
      flags$flag_show_dotplot_gmcs <- F
      if (input$I_RESULT_TABLE_GMCS_MODE_bis=="percentage"){
        for (ii in 1:length(levels(values$sample.class))) {
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii,"_bis"), min = 0, max = 100, value = c(0,100), step = 1)
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii), min = 0, max = 100, value = c(0,100), step = 1)
        }
      } else {
        for (ii in 1:length(levels(values$sample.class))) {
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii,"_bis"), min = 0,
                            max = sum(values$sample.class == levels(values$sample.class)[ii]),
                            value = c(0,sum(values$sample.class == levels(values$sample.class)[ii])), step = 1)
          updateSliderInput(session, paste0("I_perc_GEA_class_",ii), min = 0,
                            max = sum(values$sample.class == levels(values$sample.class)[ii]),
                            value = c(0,sum(values$sample.class == levels(values$sample.class)[ii])), step = 1)
        }
      }
    })
    
    
    
    # Update the table ####
    observeEvent(c(input$I_perc_GEA_class_1_bis, input$I_perc_GEA_class_2_bis,
                   input$I_perc_GEA_class_3_bis, input$I_perc_GEA_class_4_bis,
                   input$I_perc_GEA_class_5_bis, input$I_perc_GEA_class_6_bis,
                   input$I_perc_GEA_class_7_bis, input$I_perc_GEA_class_8_bis,
                   input$I_perc_GEA_class_9_bis, input$I_perc_GEA_class_10_bis,
                   input$I_perc_GEA_class_11_bis, input$I_perc_GEA_class_12_bis,
                   input$I_perc_GEA_class_13_bis, input$I_perc_GEA_class_14_bis,
                   input$I_perc_GEA_class_15_bis, input$I_perc_GEA_class_16_bis,
                   input$I_perc_GEA_class_17_bis, input$I_perc_GEA_class_18_bis,
                   input$I_perc_GEA_class_19_bis, input$I_perc_GEA_class_20_bis,
                   input$I_perc_GEA_class_21_bis, input$I_perc_GEA_class_22_bis,
                   input$I_perc_GEA_class_23_bis, input$I_perc_GEA_class_24_bis,
                   input$I_perc_GEA_class_25_bis, input$I_perc_GEA_class_26_bis,
                   input$I_perc_GEA_class_27_bis, input$I_perc_GEA_class_28_bis,
                   input$I_perc_GEA_class_29_bis, input$I_perc_GEA_class_30_bis),{
                     
                     if (length(levels(values$sample.class))>=1) values$filtersGEA[[1]] <- input$I_perc_GEA_class_1_bis
                     if (length(levels(values$sample.class))>=2) values$filtersGEA[[2]] <- input$I_perc_GEA_class_2_bis
                     if (length(levels(values$sample.class))>=3) values$filtersGEA[[3]] <- input$I_perc_GEA_class_3_bis
                     if (length(levels(values$sample.class))>=4) values$filtersGEA[[4]] <- input$I_perc_GEA_class_4_bis
                     if (length(levels(values$sample.class))>=5) values$filtersGEA[[5]] <- input$I_perc_GEA_class_5_bis
                     if (length(levels(values$sample.class))>=6) values$filtersGEA[[6]] <- input$I_perc_GEA_class_6_bis
                     if (length(levels(values$sample.class))>=7) values$filtersGEA[[7]] <- input$I_perc_GEA_class_7_bis
                     if (length(levels(values$sample.class))>=8) values$filtersGEA[[8]] <- input$I_perc_GEA_class_8_bis
                     if (length(levels(values$sample.class))>=9) values$filtersGEA[[9]] <- input$I_perc_GEA_class_9_bis
                     if (length(levels(values$sample.class))>=10) values$filtersGEA[[10]] <- input$I_perc_GEA_class_10_bis
                     if (length(levels(values$sample.class))>=11) values$filtersGEA[[11]] <- input$I_perc_GEA_class_11_bis
                     if (length(levels(values$sample.class))>=12) values$filtersGEA[[12]] <- input$I_perc_GEA_class_12_bis
                     if (length(levels(values$sample.class))>=13) values$filtersGEA[[13]] <- input$I_perc_GEA_class_13_bis
                     if (length(levels(values$sample.class))>=14) values$filtersGEA[[14]] <- input$I_perc_GEA_class_14_bis
                     if (length(levels(values$sample.class))>=15) values$filtersGEA[[15]] <- input$I_perc_GEA_class_15_bis
                     if (length(levels(values$sample.class))>=16) values$filtersGEA[[16]] <- input$I_perc_GEA_class_16_bis
                     if (length(levels(values$sample.class))>=17) values$filtersGEA[[17]] <- input$I_perc_GEA_class_17_bis
                     if (length(levels(values$sample.class))>=18) values$filtersGEA[[18]] <- input$I_perc_GEA_class_18_bis
                     if (length(levels(values$sample.class))>=19) values$filtersGEA[[19]] <- input$I_perc_GEA_class_19_bis
                     if (length(levels(values$sample.class))>=20) values$filtersGEA[[20]] <- input$I_perc_GEA_class_20_bis
                     if (length(levels(values$sample.class))>=21) values$filtersGEA[[21]] <- input$I_perc_GEA_class_21_bis
                     if (length(levels(values$sample.class))>=22) values$filtersGEA[[22]] <- input$I_perc_GEA_class_22_bis
                     if (length(levels(values$sample.class))>=23) values$filtersGEA[[23]] <- input$I_perc_GEA_class_23_bis
                     if (length(levels(values$sample.class))>=24) values$filtersGEA[[24]] <- input$I_perc_GEA_class_24_bis
                     if (length(levels(values$sample.class))>=25) values$filtersGEA[[25]] <- input$I_perc_GEA_class_25_bis
                     if (length(levels(values$sample.class))>=26) values$filtersGEA[[26]] <- input$I_perc_GEA_class_26_bis
                     if (length(levels(values$sample.class))>=27) values$filtersGEA[[27]] <- input$I_perc_GEA_class_27_bis
                     if (length(levels(values$sample.class))>=28) values$filtersGEA[[28]] <- input$I_perc_GEA_class_28_bis
                     if (length(levels(values$sample.class))>=29) values$filtersGEA[[29]] <- input$I_perc_GEA_class_29_bis
                     if (length(levels(values$sample.class))>=30) values$filtersGEA[[30]] <- input$I_perc_GEA_class_30_bis
                     
                     # flags$flag_show_heatmap_gmcs <- F
                     # flags$flag_show_dotplot_gmcs <- F
                     values$O_table_essential_gmcs_bis_rows_selected_previous_selection <- 0
                     values$O_table_essential_gmcs_rows_selected_previous_selection <- 0
                     # input$O_table_essential_gmcs_bis_rows_selected                   
                     for (ii in 1:length(levels(values$sample.class))) {
                       updateSliderInput(session, paste0("I_perc_GEA_class_",ii), value = values$filtersGEA[[ii]])
                     }
                     
                     # hide all the tables and calculate()
                     show_button_calculate_correlations()
                   })
    
    
    # Filter the data only when selected, to save computer power
    filter_DepMap_tables <- function(){
      # browser()
      
      values$DepMap.info.filtered$DepMapEssentialityByGene <- DepMap.info.all$DepMapEssentialityByGene %>%
        filter(essentiality.database==input$I_DepMap_database) %>%
        filter(ENSEMBL %in% gMCS.info$selected$genes.gMCSs.ENSEMBL)
      
      values$DepMap.info.filtered$DepMapExpressionByGene <- DepMap.info.all$DepMapExpressionByGene %>%
        filter(gMCS.database==input$I_DepMap_GMCS_LIST) %>%
        filter(UNIT==input$I_DepMap_unit) %>%
        filter(ENSEMBL %in% gMCS.info$selected$genes.gMCSs.ENSEMBL)
      
      values$DepMap.info.filtered$DepMapCorrelationByGene <- DepMap.info.all$DepMapCorrelationByGene %>%
        filter(gMCS.database==input$I_DepMap_GMCS_LIST) %>%
        filter(essentiality.database==input$I_DepMap_database) %>%
        filter(UNIT==input$I_DepMap_unit) %>%
        filter(ENSEMBL %in% gMCS.info$selected$genes.gMCSs.ENSEMBL)
      
      values$DepMap.info.filtered$DepMapGeneExpression <- DepMap.info.all$DepMapGeneExpression %>%
        filter(UNIT==input$I_DepMap_unit) %>%
        filter(ENSEMBL %in% gMCS.info$selected$genes.gMCSs.ENSEMBL)
    }
    
    # hide all the results if something has changed
    observeEvent(c(input$I_DepMap_database,
                   input$I_DepMap_database_filter_mode,
                   input$I_DepMap_database_filter_selected,
                   input$I_DepMap_unit,
                   input$I_DepMap_GMCS_LIST), {
                     if (RealTimeTables_mp5){
                       filter_DepMap_tables()
                     }
                   })
    
    # hide all the results if something has changed
    observeEvent(c(input$I_RESULT_TABLE_GMCS_MODE_bis,
                   input$I_DepMap_database,
                   input$I_DepMap_database_filter_mode,
                   input$I_DepMap_database_filter_selected,
                   input$I_DepMap_unit,
                   input$I_DepMap_GMCS_LIST,
                   input$I_DepMap_geneORgmcs), {
                     
                     if (!RealTimeTables_mp5){
                       # do nothing and show the button
                       show_button_calculate_correlations()
                       values$DepMap.info.filtered$DepMapEssentialityByGene <- NaN
                       values$DepMap.info.filtered$DepMapCorrelationByGene <- NaN
                       values$DepMap.info.filtered$DepMapGeneExpression <- NaN
                       values$DepMap.info.filtered$DepMapGeneExpression <- NaN
                     }
                   })
    
    # calculate correlation table when filters are applied ####
    ResultsTableCorrFiltersGene <- reactive({
      if (flags$flag_show_table_gmcs_correlation){
        
        # browser()
        # sizes <- sapply(reactiveValuesToList(values, all.names = T), function(n) object.size(n), simplify = FALSE); print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))
        # sizes <- sapply(gMCS.info$selected, function(n) object.size(n), simplify = FALSE); print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))
        # sizes <- sapply(reactiveValuesToList(metTasks, all.names = T), function(n) object.size(n), simplify = FALSE); print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))
        # 
        
        if (input$I_DepMap_database_filter_mode != "none" & length(setdiff(input$I_DepMap_database_filter_selected, '---'))>0) {
          
          print("Calculate new correlations based on selected data")
          
          # START <- Sys.time()
          aux1 <- unique(DepMap.info.all$dictionary.CCLE$DepMap_ID[DepMap.info.all$dictionary.CCLE[,input$I_DepMap_database_filter_mode] %in% input$I_DepMap_database_filter_selected] )
          print("Calculate new correlations based on selected data: 1/4")
          
          aux2 <- values$DepMap.info.filtered$DepMapEssentialityByGene %>%
            dplyr::filter(DepMap_ID %in% aux1)
          
          print("Calculate new correlations based on selected data: 2/4")
          # END[1] <- Sys.time()
          aux3 <- values$DepMap.info.filtered$DepMapExpressionByGene %>%
            dplyr::filter(DepMap_ID %in% aux1)
          
          print("Calculate new correlations based on selected data: 3/4")
          # END[2] <- Sys.time()
          
          # browser()
          aux <- merge(na.omit(aux2) %>% select("DepMap_ID", "essentiality_score","ENSEMBL"),
                       na.omit(aux3) %>% select("DepMap_ID", "logTPM",  "ENSEMBL")) %>%
            as_tibble() %>% 
            group_by(ENSEMBL) %>%
            # multidplyr::partition(cluster = cluster_multidplyr) %>% #cluster_library("dplyr") %>% cluster_library("rstatix") %>%
            cor_test(essentiality_score, logTPM) %>% adjust_pvalue(method = "fdr") %>%
            # collect() %>% as_tibble() %>%
            rename(ENSEMBL = ENSEMBL) %>% rename(rho = cor) %>% rename(p.value = p)  %>%
            dplyr::select(ENSEMBL, rho, p.value, p.adj)
          print("Calculate new correlations based on selected data: 4/4")
        } else {
          aux = values$DepMap.info.filtered$DepMapCorrelationByGene %>%
            dplyr::select(ENSEMBL, rho, p.value, p.adj)
        }
      } else {
        aux = NULL
      }
      return(aux)
    })
    
    
    # calculate which table to show, taking into account all gMCS at a time for the filters ####
    ResultsCorrDepMapTable <- reactive({
      # if (flags$flag_show_table_essentiality){
      if (!is.null(values$ResultsEssentiality) & flags$flag_show_table_gmcs_correlation){
        list.gMCS.essential.aux <- values$ResultsEssentiality$list.gMCS.essential
        # preprocess to obtain ratio or number
        if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
          list.gene.essential.aux <- values$ResultsEssentiality$num.essential.gene
        } else {
          list.gene.essential.aux <- values$ResultsEssentiality$ratio.essential.gene
        }
        
        if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
          for (i in 1:length(levels(values$sample.class)))
          {
            list.gMCS.essential.aux[,levels(values$sample.class)[i]] = list.gMCS.essential.aux[,levels(values$sample.class)[i]] * sum(values$sample.class == levels(values$sample.class)[i])
          }
        }
        for (i in 1:length(levels(values$sample.class)))      {
          if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
            idx <- (list.gene.essential.aux[,levels(values$sample.class)[i]] >= values$filtersGEA[[i]][1]) &
              (list.gene.essential.aux[,levels(values$sample.class)[i]] <= values$filtersGEA[[i]][2])
          } else {
            idx <- (list.gene.essential.aux[,levels(values$sample.class)[i]] >= values$filtersGEA[[i]][1]/100) &
              (list.gene.essential.aux[,levels(values$sample.class)[i]] <= values$filtersGEA[[i]][2]/100)
          }
          list.gene.essential.aux <- list.gene.essential.aux[idx,]
        }
        rownames(list.gene.essential.aux) <- 1:nrow(list.gene.essential.aux)
        
        aux <- ResultsTableCorrFiltersGene()
        
        # browser()
        
        list.gene.essential.aux <- merge(list.gene.essential.aux, aux, all.x = T)
        
        # browser()
        
        values$list.corr.essential.filtered <- list.gene.essential.aux
        if (input$I_RESULT_TABLE_GMCS_MODE=="number"){
          ResultsCorrDepMapTable <- datatable(values$list.corr.essential.filtered, filter = "top", selection = "single", rownames = F)
        } else {
          ResultsCorrDepMapTable <- datatable(values$list.corr.essential.filtered, filter = "top", selection = "single", rownames = F)   %>%  formatPercentage(columns = levels(values$sample.class),digits=2)
        }
        ResultsCorrDepMapTable <- ResultsCorrDepMapTable %>% formatRound(columns = c("rho"),digits=3) %>% formatSignif(columns = c( "p.value", "p.adj"),digits=3)
      } else {
        # browser()
        values$list.corr.essential.filtered <- merge(gMCS.info$selected$table.genes.HumanGEM, ResultsTableCorrFiltersGene())
        ResultsCorrDepMapTable <- datatable(values$list.corr.essential.filtered, filter = "top", selection = "single", rownames = F)
      }
    })
    
    
    # Show gMCS essentials
    observeEvent(c(flags$flag_show_table_gmcs_correlation, input$I_DepMap_geneORgmcs), {
      
      # browser()
      
      if (!flags$flag_show_table_gmcs_correlation){
        # print(flags$flag_show_table_essentiality)
        # print(flags$flag_show_heatmap_gmcs)
        # print(flags$flag_show_dotplot_gmcs)
        hideElement("mp4"); hideElement("mp4_heatmap_1"); hideElement("mp4_heatmap_2"); hideElement("mp4_heatmap_3")
        hideElement("mp5"); hideElement("mp5_boxplot")
        
        output$O_table_essential_genes <- NULL
        
      } else {
        showElement("mp4"); showElement("mp5")
        if (is.null(input$I_DepMap_geneORgmcs)){
          output$O_table_essential_gmcs_bis <- DT::renderDataTable({ResultsCorrDepMapTable()})
        } else if (input$I_DepMap_geneORgmcs=="gMCS") {
          output$O_table_essential_gmcs_bis <- DT::renderDataTable({ResultsGMCSsTable()})
        } else {
          output$O_table_essential_gmcs_bis <- DT::renderDataTable({ResultsCorrDepMapTable()})
        }
      }
    })
    
    
    observeEvent(input$O_table_essential_gmcs_bis_rows_selected, {
      if (values$O_table_essential_gmcs_bis_rows_selected_previous_selection!=input$O_table_essential_gmcs_bis_rows_selected) {
        
        flags$flag_show_dotplot_gmcs <- F
        print(paste0("gMCSs table row: ", input$O_table_essential_gmcs_bis_rows_selected))
        print(values$list.gMCS.essential.filtered[input$O_table_essential_gmcs_bis_rows_selected,])
        if (!is.null(values$list.gMCS.essential.filtered)){ flags$flag_show_dotplot_gmcs <- T }
        values$O_table_essential_gmcs_rows_bis_selected_previous_selection=input$O_table_essential_gmcs_bis_rows_selected
        
        flags$flag_show_heatmap_gmcs <- F
        flags$flag_show_dotplot_gmcs <- T
        
      } else {
        print(paste0("gMCSs table row: ", input$O_table_essential_gmcs_bis_rows_selected, ", same as before"))
      }
    })
    
    
    
    
    
    # generate options for filters #####
    observeEvent(input$I_DepMap_database_filter_mode, {
      
      if(input$I_DepMap_database_filter_mode == "none") {
        options_filters = c()
      } else {
        options_filters = unique(DepMap.info.all$dictionary.CCLE[,input$I_DepMap_database_filter_mode])
      }
      
      options_filters <- c("---", options_filters)
      options_filters <- as.list(options_filters)
      names(options_filters) <- unlist(options_filters)
      
      updateSelectizeInput(session, inputId = "I_DepMap_database_filter_selected",
                           choices = options_filters, selected = "---")
    })
    
    
    observeEvent(input$I_DepMap_database, {
      if(input$I_DepMap_database == "Achilles") {
        updateNumericInput(session, inputId = "I_show_DepMap_threshold", value = -0.6)
      } else {
        updateNumericInput(session, inputId = "I_show_DepMap_threshold", value = NaN)
      }
      flags$flag_show_dotplot_gmcs <- F
    })
    
    observeEvent(c(input$I_DepMap_database,
                   input$I_DepMap_database_filter_mode,
                   input$I_DepMap_unit), {
                     flags$flag_show_dotplot_gmcs <- F
                   })
    
    
    # show the selected correlation with DepMap for certain gene ####
    observeEvent(c(flags$flag_show_dotplot_gmcs,
                   input$I_DepMap_database_filter_show_only,
                   input$I_DepMap_LinearRegression,
                   input$I_show_SYMBOL_heatmap_bis,
                   input$I_show_DepMap_threshold,
                   input$I_show_separate_gMCS,
                   input$O_table_essential_gmcs_bis_rows_selected), {
                     
                     if (flags$flag_show_dotplot_gmcs){
                       # showElement("mp5_dotplot_1")
                       showElement("mp5_dotplot_2")
                       showElement("mp5_dotplot_3")
                       
                       if (is.null(input$I_DepMap_geneORgmcs)){
                         aux <- values$list.corr.essential.filtered[input$O_table_essential_gmcs_bis_rows_selected,]
                       } else if (input$I_DepMap_geneORgmcs=="gMCS") {
                         aux <- values$list.gMCS.essential.filtered[input$O_table_essential_gmcs_bis_rows_selected,]
                       } else { 
                         aux <- values$list.corr.essential.filtered[input$O_table_essential_gmcs_bis_rows_selected,]
                       }
                       print("start ShowDotplot_DepMap")
                       values$plot.obj.DepMap <- ShowDotplot_DepMap(DepMap.info.all = values$DepMap.info.filtered,
                                                                    gMCS.info.all = gMCS.info$selected,
                                                                    gmcs_database = input$I_DepMap_GMCS_LIST,
                                                                    gene.target.info = aux,
                                                                    database = input$I_DepMap_database,
                                                                    database_filter_mode = input$I_DepMap_database_filter_mode,
                                                                    database_filter_selected = input$I_DepMap_database_filter_selected,
                                                                    flag_database_filter_show_only = ifelse(input$I_DepMap_database_filter_mode=="none", 0,input$I_DepMap_database_filter_show_only),
                                                                    flag_LinearRegression = input$I_DepMap_LinearRegression,
                                                                    flag_show_SYMBOL = input$I_show_SYMBOL_heatmap_bis,
                                                                    threshold_value = input$I_show_DepMap_threshold,
                                                                    flag_color_by_gMCS = input$I_show_separate_gMCS,
                                                                    database_unit = input$I_DepMap_unit,
                                                                    by_gMCS = input$I_DepMap_geneORgmcs=="gMCS")
                       
                       output$O_dotplot_seleted_gmcs <- renderPlot({ggarrange(values$plot.obj.DepMap)})
                       
                       print("end ShowDotplot_DepMap")
                     } else {
                       # hideElement("mp5_dotplot_1")
                       hideElement("mp5_dotplot_2")
                       hideElement("mp5_dotplot_3")
                       
                       values$plot.obj.DepMap <- ggplot() + theme_void()
                       output$O_dotplot_seleted_gmcs <- renderPlot({ggplot() + theme_void()})
                     }
                   })
    
    
    # Save the results ####
    output$O_dotplot_seleted_gmcs_download_1 <- downloadHandler(filename = paste0("gMCStool_dotplot_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".png"),
                                                                content = function(fname){ggsave(values$plot.obj.DepMap, file = fname,
                                                                                                 width = 10, height = 8,
                                                                                                 units = "in", dpi = 400, bg = "white")})
    output$O_dotplot_seleted_gmcs_download_2 <- downloadHandler(filename = paste0("gMCStool_dotplot_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".jpg"),
                                                                content = function(fname){ggsave(values$plot.obj.DepMap, file = fname,
                                                                                                 width = 10, height = 8,
                                                                                                 units = "in", dpi = 400, bg = "white")})
    output$O_dotplot_seleted_gmcs_download_3 <- downloadHandler(filename = paste0("gMCStool_dotplot_", format(Sys.time(), "%Y_%m_%d_%Hh%Mm"), ".RDS"),
                                                                content = function(fname){saveRDS(values$plot.obj.DepMap, file = fname)})
    
    
  }
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # END OF FILE ####
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  output$O_help_download <- downloadHandler(filename = "gMCStool_help_tab.pdf",
                                            content = function(fname){file.copy(from = "instructions/gMCStool_help_tab.pdf", to = fname, overwrite = T)})
  
  output$iframe <- renderUI({
    tags$iframe(src='https://www.google.com/maps/embed?pb=!1m14!1m8!1m3!1d11616.49197729335!2d-1.9834123!3d43.2907459!3m2!1i1024!2i768!4f13.1!3m3!1m2!1s0x0%3A0x1da2ea8909de8dab!2sCEIT-Tecnun.%20Sede%20Miram%C3%B3n%20UNAV!5e0!3m2!1ses!2ses!4v1572350894436!5m2!1ses!2ses', width='70%', height='500px', frameborder='0', style='border:0;', allowfullscreen=FALSE)
    # tags$iframe(src='https://www.google.com/maps/place/CEIT-Tecnun.+Sede+Miram%C3%B3n+UNAV/@43.2907654,-1.9837456,17.87z/data=!4m5!3m4!1s0xd51afde7e1421cd:0x1da2ea8909de8dab!8m2!3d43.2907482!4d-1.9834039', width='70%', height='500px', frameborder='0', style='border:0;', allowfullscreen=FALSE)
  })
  
}



# Run the application 
app <- shinyApp(ui = ui, server = server)

