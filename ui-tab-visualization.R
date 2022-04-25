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



tabPanel(title = "4. Visualization",
         
         # useShinyjs(),  # Set up shinyjs
         
         sidebarLayout(
           # Sidebar
           sidebarPanel(           
             width = 3, # 3 Ferca value"
             h4("Your Analysis"),
             h5("Selected gMCS database:"),
             h6(textOutput("O_txtout_I_GMCS_LIST")),
             h5("Selected thresholding method:"),
             h6(textOutput("O_txtout_I_TH_METHOD")),
             #br(),
             hr(style="border-color: black; background-color: black;color: black;"),            
            
             
             h4("Define the target classes:"),
             checkboxGroupInput("I_isTarget_GEA_class", 
                                NULL, #"Check only the desired target classes:",
                                c("val1","val2"),
                                selected = c("val1","val2")),
             
             
             h4("Percentage of gene essentiality in target vs non-target class"),
             sliderInput(inputId = "I_perc_GEA_class_target_filter", h5("Threshold values for essentiality in target class (disease)"), min = 0, max = 100, value = c(1,100), dragRange = F, step = 1),
             sliderInput(inputId = "I_perc_GEA_class_non_target_filter", h5("Threshold values for essentiality in non-target class (healthy)"), min = 0, max = 100, value = c(0,0), dragRange = F, step = 1),
             br(),
             
             fluidRow(
               column(width = 5, actionButton("I_class_1_30_filter_targeted_class",
                                              "Update all",
                                              width = '100%',heigth = '100%',
                                              style="color: #ffffff; background-color: #333333; border-color: #ffffff")),
               column(width = 2),
               column(width = 5, actionButton("I_class_1_30_unfilter_targeted_class",
                                              "Reset all",
                                              width = '100%',heigth = '100%',
                                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"))
             ),
             
             br(),
             
             
             # Slider for Gene Essentiality Analysis filtering at gene level
             h4("Gene Essentiality Filtering"),
             h5("Percentage of gene essentiality across samples in class: "),
             div(id = "mp4_class_1", sliderInput(inputId = "I_perc_GEA_class_1", textOutput("O_name_GEA_class_1"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_2", sliderInput(inputId = "I_perc_GEA_class_2", textOutput("O_name_GEA_class_2"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_3", sliderInput(inputId = "I_perc_GEA_class_3", textOutput("O_name_GEA_class_3"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_4", sliderInput(inputId = "I_perc_GEA_class_4", textOutput("O_name_GEA_class_4"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_5", sliderInput(inputId = "I_perc_GEA_class_5", textOutput("O_name_GEA_class_5"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_6", sliderInput(inputId = "I_perc_GEA_class_6", textOutput("O_name_GEA_class_6"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_7", sliderInput(inputId = "I_perc_GEA_class_7", textOutput("O_name_GEA_class_7"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_8", sliderInput(inputId = "I_perc_GEA_class_8", textOutput("O_name_GEA_class_8"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_9", sliderInput(inputId = "I_perc_GEA_class_9", textOutput("O_name_GEA_class_9"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_10", sliderInput(inputId = "I_perc_GEA_class_10", textOutput("O_name_GEA_class_10"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_11", sliderInput(inputId = "I_perc_GEA_class_11", textOutput("O_name_GEA_class_11"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_12", sliderInput(inputId = "I_perc_GEA_class_12", textOutput("O_name_GEA_class_12"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_13", sliderInput(inputId = "I_perc_GEA_class_13", textOutput("O_name_GEA_class_13"),min = 0, max = 130, value = c(0,130), dragRange = F, step = 1)),  
             div(id = "mp4_class_14", sliderInput(inputId = "I_perc_GEA_class_14", textOutput("O_name_GEA_class_14"),min = 0, max = 140, value = c(0,140), dragRange = F, step = 1)),  
             div(id = "mp4_class_15", sliderInput(inputId = "I_perc_GEA_class_15", textOutput("O_name_GEA_class_15"),min = 0, max = 150, value = c(0,150), dragRange = F, step = 1)),  
             div(id = "mp4_class_16", sliderInput(inputId = "I_perc_GEA_class_16", textOutput("O_name_GEA_class_16"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_17", sliderInput(inputId = "I_perc_GEA_class_17", textOutput("O_name_GEA_class_17"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_18", sliderInput(inputId = "I_perc_GEA_class_18", textOutput("O_name_GEA_class_18"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_19", sliderInput(inputId = "I_perc_GEA_class_19", textOutput("O_name_GEA_class_19"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_20", sliderInput(inputId = "I_perc_GEA_class_20", textOutput("O_name_GEA_class_20"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_21", sliderInput(inputId = "I_perc_GEA_class_21", textOutput("O_name_GEA_class_21"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_22", sliderInput(inputId = "I_perc_GEA_class_22", textOutput("O_name_GEA_class_22"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_23", sliderInput(inputId = "I_perc_GEA_class_23", textOutput("O_name_GEA_class_23"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_24", sliderInput(inputId = "I_perc_GEA_class_24", textOutput("O_name_GEA_class_24"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_25", sliderInput(inputId = "I_perc_GEA_class_25", textOutput("O_name_GEA_class_25"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_26", sliderInput(inputId = "I_perc_GEA_class_26", textOutput("O_name_GEA_class_26"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_27", sliderInput(inputId = "I_perc_GEA_class_27", textOutput("O_name_GEA_class_27"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_28", sliderInput(inputId = "I_perc_GEA_class_28", textOutput("O_name_GEA_class_28"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_29", sliderInput(inputId = "I_perc_GEA_class_29", textOutput("O_name_GEA_class_29"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
             div(id = "mp4_class_30", sliderInput(inputId = "I_perc_GEA_class_30", textOutput("O_name_GEA_class_30"),min = 0, max = 100, value = c(0,100), dragRange = F, step = 1)),  
           ),
           
           
           
           
           
           # Main panel where all the info/results are shown ------------------------------------------------------------------------------------------------------
           mainPanel(width = 9,
                     # tabsetPanel(
                     h3("Results of Gene Essentiality Analysis:"),
                     br(),
                     radioButtons(inputId = "I_RESULT_TABLE_GMCS_MODE", label = "Select:", 
                                  choiceValues = list("number",
                                                      "percentage"), 
                                  choiceNames =  list("number",
                                                      "percentage"), 
                                  selected = "percentage", inline = TRUE),
                     br(),
                     div(id = "mp4_button",
                         br(), br(), br(), br(),
                         fluidRow(column(width = 3),
                                  column(width = 6,
                                         actionButton("I_mp4_calculate", 
                                                      label = h4("Show result table"), width = '100%',heigth = '200%',
                                                      style=STYLE_NEXT_BUTTONS), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                                  ),
                                  column(width = 3)),
                         br(), br(), br(), br(),
                     ),
                     
                     div(id = "mp4",
                         div(withSpinner(DT::dataTableOutput(outputId = "O_table_essential_gmcs")), style = "font-size: 95%; width: 100%"),
                         br(),
                     ),
                     div(id = "mp4_heatmap",
                         br(),
                         fluidRow(
                           column(width = 3, checkboxInput(inputId = "I_include_boxplots_heatmap",
                                                           label = "Include boxplots of expression in heatmap?", value = TRUE)),
                           column(width = 3, checkboxInput(inputId = "I_show_SYMBOL_heatmap",
                                                           label = "Show gene SYMBOL in heatmap? (Default is ENSEMBL)", value = TRUE)),
                           column(width = 3, selectInput("I_colors_heatmap", h5("Select colors pallete for the heatmap"),
                                                         choices = list("Blue -- white -- red" = "royalblue1_white_firebrick1",
                                                                        "Purple -- white -- orange" = "#902EBB_white_#F4831B"),
                                                         selected = "royalblue1_white_firebrick1")),
                           column(width = 3, selectInput("I_colors_colors_annotation", h5("Select colors pallete for the annotation"),
                                                         choices = list("rainbow" = 1,
                                                                        "hue pal" = 2),
                                                         selected = 1))
                         ),
                         br(),
                     ),
                     div(id = "mp4_heatmap_1",
                         h5("Target metabolic tasks:"),
                         div(tableOutput(outputId = "O_table_essential_gmcs_selected_tasks"), style = "font-size: 95%; width: 100%"),
                         h5("Target metabolite (only active for biomass production):"),
                         div(tableOutput(outputId = "O_table_essential_gmcs_selected_mets"), style = "font-size: 95%; width: 100%"),
                         br(),
                         plotOutput(outputId = "O_heatmap_seleted_gmcs", height="600px", width="1200px"),
                         br()
                     ),
                     div(id = "mp4_heatmap_2",
                         plotOutput(outputId = "O_boxplot_seleted_gmcs", height="400px", width="1200px"),
                         br()
                     ),
                     div(id = "mp4_heatmap_3",
                         fluidRow(
                           column(width = 4, downloadButton("O_heatmap_seleted_gmcs_download_1","Download plot (PNG)")),
                           column(width = 4, downloadButton("O_heatmap_seleted_gmcs_download_2","Download plot (JPG)")),
                           column(width = 4, downloadButton("O_heatmap_seleted_gmcs_download_3","Download ggplot (RDS)")),
                         )
                     ),
                     
                     div(id = "mp4_loading",
                         br(),
                         fluidRow(column(width = 3),
                                  column(width = 6, img(src="loading-74.gif", align = "center",height='350px',width='350px')),
                                  column(width = 3)),
                         br()
                     )
           )   
         )
)
