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


tabPanel(title = "3. Predict Essential Genes and DKOs",
         
         sidebarLayout(
           # Sidebar
           sidebarPanel(width = 3,
                        
                        # Select list of gMCSs to use in the study
                        
                        h3("Select the gMCS database: "),
                        
                        radioButtons(inputId = "I_GMCS_1_EssentialTask_or_Only", label = 'Select one:' , 
                                     choiceValues = list("EssentialTasks",
                                                         "Only",
                                                         "Custom"), 
                                     choiceNames =  list("All metabolic essential tasks",
                                                         "Only biomass production",
                                                         "Selected tasks"), 
                                     selected = "EssentialTasks"),
                        
                        radioButtons(inputId = "I_GMCS_1_CultureMedia_or_Full", label = 'Select biomass production restrictions:', 
                                     choiceValues = list("CultureMedium",
                                                         "FullMedium"), 
                                     choiceNames =  list("Growth on Ham's medium",
                                                         "Growth on unconstrained medium (all uptakes available)"), 
                                     selected = "CultureMedium"),
                        
                        
                        # Gene expression thresholding methodology
                        
                        h3("Gene expression thresholding method"),
                        radioButtons(inputId = "I_TH_METHOD", label = "Select one:", 
                                     # change the comments to implement singleTH
                                     # choiceValues = list("gmcsTH", "localT2","singleTH", "fastcormics"),
                                     # choiceNames =  list("gmcsTH", "localT2","singleTH", "fastcormics"),
                                     choiceValues = list("gmcsTH", "localT2"),
                                     choiceNames =  list("gmcsTH", "localT2"),
                                     selected = "gmcsTH"),
                        
                        # depending on the thresholding technology, show a different menu
                        div(id = "I_perc_gmcsTH_show",
                            fluidRow(
                              column(width = 9, sliderInput(inputId = "I_perc_gmcsTH", "percentile [%] of expression threshold",
                                                            min = 0, max = 1, value = c(0.05), step = 0.01)),
                              column(width = 3, br(), br(), textInput(inputId = "I_perc_gmcsTH_txt", NULL, value = "0.05"))
                            )
                        ),
                        div(id = "I_perc_localT2_show", 
                            radioButtons(inputId = "I_perc_localT2", label = "Select:", 
                                         choiceValues = list("all_genes_HumanGEM",
                                                             "all_genes_gMCSs"), 
                                         choiceNames =  list("all genes in HumanGEM",
                                                             "all genes in selected gMCSs database"), 
                                         selected = "all_genes_gMCSs")#,
                            # h5("all samples are considered to be from the same cohort")
                        ),
                        div(id = "I_perc_singleTH_show", 
                            numericInput(inputId = "I_perc_singleTH", 
                                         label = "Numeric value of expression threshold", 
                                         value = 1)
                        ),
                        div(id = "I_perc_fastcormics_show", 
                            checkboxInput(inputId = "I_perc_fastcormics", 
                                          label = "Should the data be log2 transformed for the threshold?", 
                                          value = TRUE)
                        ),
                        
                        
                        div(id = "I_hide_advance_settings_calculate_GEA", 
                            actionButton("I_action_show_advance_settings_calculate_GEA",
                                         "Advance settings",
                                         width = '100%',heigth = '100%',
                                         style="color: #ffffff; background-color: #333333; border-color: #ffffff")
                        ),
                        
                        div(id = "I_show_advance_settings_calculate_GEA", 
                            actionButton("I_action_hide_advance_settings_calculate_GEA",
                                         "Hide advance settings",
                                         width = '100%',heigth = '100%',
                                         style="color: #ffffff; background-color: #333333; border-color: #ffffff"),
                            
                            selectInput("I_GEA_settings_simplify", 
                                        h5("Select the simplification of gMCS along tasks:"),
                                        choices = list("None simplification" = 0,
                                                       "Simple simplification" = 1,
                                                       "Complete simplification" = 2),
                                        selected = "None simplification"),
                            
                            numericInput("I_GEA_settings_threshold_logFC", 
                                         h5("Select threshold_logFC to consider a gene as ON:"),
                                         1e-3),
                            
                            numericInput("I_GEA_settings_nWorkers", 
                                         h5("Select the number of cores:"),
                                         nWorkers, min = 0, max = nWorkers),
                            
                        ),
                        
                        
                        br(),
                        actionButton("I_action_predictEssentialGene", "Calculate single KO!", width = '100%',heigth = '200%',
                                     style=STYLE_NEXT_BUTTONS), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                        tags$head(tags$style(type = "text/css",".btn-default { background-image: none;}")),
                        tags$head(tags$style(type = "text/css",".btn-default:hover { background-image: none;}")),
                        
                        br(), br(), 
                        actionButton("I_action_predictEssentialPairGene", "Calculate double KO!", width = '100%',heigth = '200%',
                                     style=STYLE_NEXT_BUTTONS), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                        tags$head(tags$style(type = "text/css",".btn-default { background-image: none;}")),
                        tags$head(tags$style(type = "text/css",".btn-default:hover { background-image: none;}")),
                        
           ),
           # Main panel
           mainPanel(# co_FDR = 0.4, co_log2FC = 2
             # tabsetPanel(
             # h3("Ranking of Essential Genes",
             
             h3("Load examples or previously calculated results:"),
             
             fluidRow(
               column(width = 6, actionButton("I_ResultsEssentiality_show_input_RDATA",
                                              "Load previously calculated results",
                                              width = '100%',heigth = '100%',
                                              style="color: #ffffff; background-color: #888888; border-color: #555555")),
               column(width = 6, actionButton("I_ResultsEssentiality_show_input_examples",
                                              "Load Example Results",
                                              width = '100%',heigth = '100%',
                                              style="color: #ffffff; background-color: #888888; border-color: #555555")),
             ),
             
             div(id = "mp3_ResultsEssentiality_input_RDATA", 
                 br(),
                 fileInput("ResultsEssentiality_input", "Choose .RData file to upload precomputed results",
                           width = "100%",
                           accept = c(
                             ".RData",
                             ".Rdata")
                 )
             ),
             
             div(id = "mp3_ResultsEssentiality_input_example", 
                 h4("Select example to upload"),
                 actionButton("I_ResultsEssentiality_input_example_gmcsTH_TPM",
                              "Load Example Data precomputed Results for gmcsTH(5%) [TPM]",
                              width = '100%',heigth = '100%',
                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                 actionButton("I_ResultsEssentiality_input_example_gmcsTH_log2TPM",
                              "Load Example Data precomputed Results for gmcsTH(5%) [log2(TPM+1)]",
                              width = '100%',heigth = '100%',
                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                 
                 actionButton("I_ResultsEssentiality_input_example_localT2_TPM",
                              "Load Example Data precomputed Results for localT2 (all genes in gMCS) [TPM]",
                              width = '100%',heigth = '100%',
                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                 actionButton("I_ResultsEssentiality_input_example_localT2_log2TPM",
                              "Load Example Data precomputed Results for localT2 (all genes in gMCS) [log2(TPM+1)]",
                              width = '100%',heigth = '100%',
                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                 
                 actionButton("I_ResultsEssentiality_input_all_DepMap_gmcsTH_TPM",
                              "Load Results of essentiality analysis for all cell lines in DepMap, gmcsTH(5%) [TPM]",
                              width = '100%',heigth = '100%',
                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
             ),
             br(),
             div(id = "mp3",
                 h3("Results of Gene Essentiality Analysis:"),
                 br(),
                 fluidRow(
                   column(width = 6, selectInput(inputId = "I_RESULT_TABLE_SINGLE_DOUBLE",
                                                 label = "Select which results to show:", 
                                                 choices =  list("Single KO", "Double KO (only)", "Double KO + single KOs"), 
                                                 selected = "Single KO")
                   ),
                   column(width = 6, radioButtons(inputId = "I_RESULT_TABLE_MODE", 
                                                  label = "Select visualization mode:", 
                                                  choiceValues = list("number",
                                                                      "percentage"), 
                                                  choiceNames =  list("number",
                                                                      "percentage"), 
                                                  selected = "percentage", inline = TRUE)
                   ),
                 ),
                 
                 br(),
                 div(withSpinner(DT::dataTableOutput(outputId = "O_table_essential_genes")), style = "font-size: 95%; width: 100%"),
                 br(),
                 fluidRow(column(width = 4, downloadButton("O_results_essential_genes_simple_download","Download table")),
                          column(width = 4, downloadButton("O_results_essential_genes_all_download","Download Excel (all)")),
                          column(width = 4, downloadButton("O_results_all_download","Download Rdata"))
                 ),
                 br(),
                 actionButton("I_mp3_next", "Next!", width = '100%',heigth = '200%',
                              style=STYLE_NEXT_BUTTONS), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                 br(),
             ),
             div(id = "mp3_loading",
                 h3("Performing the Gene Essentiality Analysis:"),
                 br(),
                 fluidRow(column(width = 3),
                          column(width = 6, img(src="loading-74.gif", align = "center",height='350px',width='350px')),
                          column(width = 3)),
                 br()
             )
             
           )
         )
)


