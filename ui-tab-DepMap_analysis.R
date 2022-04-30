## ==================================================================================== ##
# gMCStool Shiny App for predicting of essential genes.
# Copyright (C) 2021 Luis V. Valcarcel and Francisco J. Planes
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




tabPanel(title = "5. DepMap Analysis",
         
         # useShinyjs(),  # Set up shinyjs
         
         sidebarLayout(
           # Sidebar
           sidebarPanel(           
             width = 3, # 3 Ferca value"
             h4("Your Analysis"),
             h5("Selected gMCS database:"),
             h6(textOutput("O_txtout_I_GMCS_LIST_bis")),
             h5("Selected thresholding method:"),
             h6(textOutput("O_txtout_I_TH_METHOD_bis")),
             #br(),
             hr(style="border-color: black; background-color: black;color: black;"),
             # br(),
             
             # Slider for Gene Essentiality Analysis filtering at gene level
             # h4("Gene Essentiality Analysis filtering at gene level (ratio)"),
             h4("Gene Essentiality Filtering"),
             h5(textOutput("O_txtout_title_filter_bis")),
             # h5("Percentage of gene essentiality across samples in class: "),
             div(id = "mp5_class_1", sliderInput(inputId = "I_perc_GEA_class_1_bis", textOutput("O_name_GEA_class_1_bis"), min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_2", sliderInput(inputId = "I_perc_GEA_class_2_bis", textOutput("O_name_GEA_class_2_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_3", sliderInput(inputId = "I_perc_GEA_class_3_bis", textOutput("O_name_GEA_class_3_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_4", sliderInput(inputId = "I_perc_GEA_class_4_bis", textOutput("O_name_GEA_class_4_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_5", sliderInput(inputId = "I_perc_GEA_class_5_bis", textOutput("O_name_GEA_class_5_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_6", sliderInput(inputId = "I_perc_GEA_class_6_bis", textOutput("O_name_GEA_class_6_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_7", sliderInput(inputId = "I_perc_GEA_class_7_bis", textOutput("O_name_GEA_class_7_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_8", sliderInput(inputId = "I_perc_GEA_class_8_bis", textOutput("O_name_GEA_class_8_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_9", sliderInput(inputId = "I_perc_GEA_class_9_bis", textOutput("O_name_GEA_class_9_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_10", sliderInput(inputId = "I_perc_GEA_class_10_bis", textOutput("O_name_GEA_class_10_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_11", sliderInput(inputId = "I_perc_GEA_class_11_bis", textOutput("O_name_GEA_class_11_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_12", sliderInput(inputId = "I_perc_GEA_class_12_bis", textOutput("O_name_GEA_class_12_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_13", sliderInput(inputId = "I_perc_GEA_class_13_bis", textOutput("O_name_GEA_class_13_bis"),min = 0, max = 130, value = c(0,130), dragRange = F, step = 1)),  
             div(id = "mp5_class_14", sliderInput(inputId = "I_perc_GEA_class_14_bis", textOutput("O_name_GEA_class_14_bis"),min = 0, max = 140, value = c(0,140), dragRange = F, step = 1)),  
             div(id = "mp5_class_15", sliderInput(inputId = "I_perc_GEA_class_15_bis", textOutput("O_name_GEA_class_15_bis"),min = 0, max = 150, value = c(0,150), dragRange = F, step = 1)),  
             div(id = "mp5_class_16", sliderInput(inputId = "I_perc_GEA_class_16_bis", textOutput("O_name_GEA_class_16_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_17", sliderInput(inputId = "I_perc_GEA_class_17_bis", textOutput("O_name_GEA_class_17_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_18", sliderInput(inputId = "I_perc_GEA_class_18_bis", textOutput("O_name_GEA_class_18_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_19", sliderInput(inputId = "I_perc_GEA_class_19_bis", textOutput("O_name_GEA_class_19_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_20", sliderInput(inputId = "I_perc_GEA_class_20_bis", textOutput("O_name_GEA_class_20_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_21", sliderInput(inputId = "I_perc_GEA_class_21_bis", textOutput("O_name_GEA_class_21_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_22", sliderInput(inputId = "I_perc_GEA_class_22_bis", textOutput("O_name_GEA_class_22_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_23", sliderInput(inputId = "I_perc_GEA_class_23_bis", textOutput("O_name_GEA_class_23_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_24", sliderInput(inputId = "I_perc_GEA_class_24_bis", textOutput("O_name_GEA_class_24_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_25", sliderInput(inputId = "I_perc_GEA_class_25_bis", textOutput("O_name_GEA_class_25_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_26", sliderInput(inputId = "I_perc_GEA_class_26_bis", textOutput("O_name_GEA_class_26_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_27", sliderInput(inputId = "I_perc_GEA_class_27_bis", textOutput("O_name_GEA_class_27_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_28", sliderInput(inputId = "I_perc_GEA_class_28_bis", textOutput("O_name_GEA_class_28_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_29", sliderInput(inputId = "I_perc_GEA_class_29_bis", textOutput("O_name_GEA_class_29_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp5_class_30", sliderInput(inputId = "I_perc_GEA_class_30_bis", textOutput("O_name_GEA_class_30_bis"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
           ),
           
           
           # Main panel where all the info/results are shown ------------------------------------------------------------------------------------------------------
           mainPanel(width = 9,
                     # tabsetPanel(
                     h3("Results of Gene Essentiality Analysis:"),
                     br(),
                     fluidRow(
                       column(width = 3,radioButtons(inputId = "I_RESULT_TABLE_GMCS_MODE_bis", h5("Select numeric display:"), 
                                                     choiceValues = list("number",
                                                                         "percentage"), 
                                                     choiceNames =  list("number",
                                                                         "percentage"), 
                                                     selected = "percentage", inline = TRUE)),
                       column(width = 3, selectInput("I_DepMap_unit", 
                                                     h5("Select the RNA-seq unit:"),
                                                     choices = list("TPM [log2(TPM+1)]" = "log2(TPM+1)",
                                                                    # "pTPM [log2(pTPM+1)]" = "log2(pTPM+1)",
                                                                    "TPM (z-scores)" = "zscores(TPM+1))",
                                                                    "log2TPM (z-scores)" = "zscores(log2(TPM+1))"),
                                                     
                                                     selected = "TPM")),
                       
                       column(width = 4, selectInput("I_DepMap_GMCS_LIST", 
                                                     h5("Select gMCS dataset:"),
                                                     choices =  list("All metabolic essential tasks with Growth on Ham's medium" = "EssentialTasks_CultureMedium",
                                                                         "All metabolic essential tasks with Growth on unconstrained medium (all uptakes available)" = "EssentialTasks_FullBiomassMedium",
                                                                         "Only biomass production on Ham's medium" = "Only_CultureMedium",
                                                                         "Only biomass production on unconstrained medium (all uptakes available)" = "Only_FullMedium"),
                                                     selected = "EssentialTasks_CultureMedium",
                                                     width = "100%")),
                       
                        # Removed due to excesive computation power
                        column(width = 2,radioButtons(inputId = "I_DepMap_geneORgmcs",
                                                      h5("Show correlation by single gMCS or aggregated:"),
                                                      choiceValues = list("gene", "gMCS"),
                                                      choiceNames =  list("aggregated", "single"),
                                                      selected = "gene", inline = TRUE)),
                     ),
                     br(),
                     
                     fluidRow(
                       column(width = 4, selectInput("I_DepMap_database", h5("Select the essentiality database you want to select"),
                                                     choices = list("CRISPR (Achilles - CERES)" = "Achilles",
                                                                    "CRISPR (Achilles - Chronos)" = "AchillesChronos",
                                                                    "CRISPR (Achilles + Sanger SCORE)" = "CRISPR",
                                                                    "DEMETER v2" = "Demeter"),
                                                     selected = "Achilles")),
                       column(width = 4, selectInput("I_DepMap_database_filter_mode", 
                                                     h5("Select the filter"),
                                                     choices = list("none" = "none",
                                                                    "Primary disease" = "primary_disease",
                                                                    "Subtype" = "Subtype",
                                                                    "Cell lineage" = "lineage"),
                                                     selected = "none")),
                       column(width = 4, selectizeInput("I_DepMap_database_filter_selected",
                                                        h5("Select the category within the filter"),
                                                        choices = list("---" = "---"),
                                                        selected = "---", multiple = T))
                     ),
                     br(),
                     
                     div(id = "mp5_button",
                         br(), br(), br(), br(),
                         
                         fluidRow(column(width = 1),
                                  column(width = 6,
                                         actionButton("I_mp5_calculate", 
                                                      label = h4("Show result table"), width = '100%',heigth = '200%',
                                                      style=STYLE_NEXT_BUTTONS), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                                         ),
                                  column(width = 3)),
                         br(), br(), br(), br(),
                     ),
                     
                     
                     div(id = "mp5",
                         div(withSpinner(DT::dataTableOutput(outputId = "O_table_essential_gmcs_bis")), style = "font-size: 95%; width: 100%"),
                         br(),
                         h4("Visualization options:"),
                         fluidRow(
                           column(width = 2, checkboxInput(inputId = "I_DepMap_database_filter_show_only", 
                                                           label = "Only show selected in filter", 
                                                           value = FALSE)),
                           column(width = 2, checkboxInput(inputId = "I_DepMap_LinearRegression", 
                                                           label = "Linear Regression?", 
                                                           value = TRUE)),
                           column(width = 3, checkboxInput(inputId = "I_show_SYMBOL_heatmap_bis",
                                                           label = "Show gene SYMBOL? (Default is ENSEMBL)",
                                                           value = T)),
                           column(width = 2, checkboxInput(inputId = "I_show_separate_gMCS",
                                                           label = "Color by gMCS?",
                                                           value = FALSE)),
                           column(width = 3, numericInput(inputId = "I_show_DepMap_threshold",
                                                          label = "Show line for essentiality threshold?",
                                                          value = NULL))
                         ),
                         br(),
                     ),
                     div(id = "mp5_dotplot_2",
                         br(),
                         plotOutput(outputId = "O_dotplot_seleted_gmcs", height="600px", width="1200px"),
                         br()
                     ),
                     div(id = "mp5_dotplot_3",
                         br(),
                         fluidRow(
                           column(width = 4, downloadButton("O_dotplot_seleted_gmcs_download_1","Download plot (PNG)")),
                           column(width = 4, downloadButton("O_dotplot_seleted_gmcs_download_2","Download plot (JPG)")),
                           column(width = 4, downloadButton("O_dotplot_seleted_gmcs_download_3","Download ggplot (RDS)")),
                         )
                     )
           )
         )   
)

