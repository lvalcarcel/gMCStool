


# CRAN's packages
pkgsNeeded <- c('shiny','shinyjs','shinythemes','shinyWidgets', 'rmarkdown',
                'shinycssloaders', 'shinydashboard', 'shinydashboardPlus',
                'shinyalert', 'shinyhelper', 'DT', 'shinyTree', 
                
                'foreach', 'doParallel', 'parallel', 'doSNOW', 'doParallel', 'snow', 'pbapply',
                'data.table', 'dplyr', 'rstatix',
                'Matrix', 'reshape2', 'gridExtra', 'pheatmap', 
                'ggplot2', 'ggpubr', 'ggrepel', 
                'scales', 'colorspace', 'circlize', 
                # 'BiocManager', 
                'openxlsx', 'viridis', 'latex2exp', 'benchmarkme',
                
                'tidyverse', 'gRbase', 'Rfast2', 'parallel', 'ParallelLogger',
                'pbmcapply',
                'doParallel', 'foreach', 'parallelly', 'future', 'future.apply') 
                

pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages())
if (length(pkgsNeeded[!pkgsInstalled]) > 0){
  install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = "https://cloud.r-project.org")
}

# Bioconductor's packages
pkgsNeeded <- c("gRbase", "graph")
pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages())
if (length(pkgsNeeded[!pkgsInstalled]) > 0){
  BiocManager::install(pkgsNeeded[!pkgsInstalled])
}


library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyalert)
library(shinycssloaders)
library(shinyTree)

library(DT)

library(foreach)
library(doParallel)
library(parallel)
library(data.table)

library(Matrix)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(openxlsx)
library(doSNOW)
library(snow)
library(viridis)
library(rstatix)
library(dplyr)
#library(feather)

library(tidyverse)
library(gRbase)
library(Rfast2)
library(parallel)
library(ParallelLogger)
library(doParallel)
library(foreach)
library(parallelly)
library(future)
library(future.apply)





# load all the scrips needed in the server #### 
# source("fun-CalculateEssentialGenes_gmcsTH.R")
# source("fun-CalculateEssentialGenes_localT2.R")
# source("fun-CalculateEssentialGenes_singleTH.R")

source("fun-CalculateEssentialGenes.R")
source("fun-CalculateGenesOnOff.R")
source("fun-SimplifyGMCS.R")

source("fun-SaveResultsInExcel.R")
source("fun-ShowHeatmap.R")


source("fun-ShowDotplot_DepMap.R")






# cat("\f")
