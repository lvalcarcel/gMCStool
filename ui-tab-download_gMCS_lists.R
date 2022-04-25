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


tabPanel(title = "1. gMCS database",
         
         sidebarLayout(
           # Sidebar
           sidebarPanel(width = 3,
                        
                        # Select list of gMCSs to use in the study
                        
                        h3("Select the  gMCS database: "),
                        # prettyRadioButtons(inputId = "I_gMCS_GMCS_LIST", label = ' ', 
                        #              choiceValues = list("EssentialTasks_CultureMedium",
                        #                                  "EssentialTasks_FullBiomass",
                        #                                  "Only_CultureBiomass",          
                        #                                  "Only_FullBiomass"  ), 
                        #              choiceNames =  list("Essential Tasks and proliferation with Culture media",
                        #                                  "Essential Tasks and proliferation with Full media",
                        #                                  "Only proliferation with Culture media",          
                        #                                  "Only proliferation with Full media"  ), 
                        #              selected = "EssentialTasks_CultureMedium",
                        #              shape = "square" ),
                        
                        # fluidRow(column(width = 6,
                        radioButtons(inputId = "I_GMCS_2_EssentialTask_or_Only", label = 'Select parameters:' , 
                                     choiceValues = list("EssentialTasks",
                                                         "Only",
                                                         "Custom"), 
                                     choiceNames =  list("All metabolic essential tasks",
                                                         "Only biomass production",
                                                         "Selected tasks"), 
                                     selected = "EssentialTasks"),
                        # ),
                        # column(width = 6, 
                        radioButtons(inputId = "I_GMCS_2_CultureMedia_or_Full", label = 'Select biomass production restrictions:', 
                                     choiceValues = list("CultureMedium",
                                                         "FullMedium"), 
                                     choiceNames =  list("Growth on Ham's medium",
                                                         "Growth on unconstrained medium (all uptakes available)"), 
                                     selected = "CultureMedium"),
                        # )
                        # ),
                        br(),
                        h4("Select the metabolic tasks within the gMCS database: "),
                        
                        fluidRow(
                          column(width = 5, actionButton("I_gMCS_GMCS_LIST_TASKS_all",
                                                         "Select all",
                                                         width = '100%',heigth = '100%',
                                                         style="color: #ffffff; background-color: #333333; border-color: #ffffff")),
                          column(width = 2),
                          column(width = 5, actionButton("I_gMCS_GMCS_LIST_TASKS_none",
                                                         "Deselect all",
                                                         width = '100%',heigth = '100%',
                                                         style="color: #ffffff; background-color: #333333; border-color: #ffffff"))
                        ),
                        
                        
                        # div(id = "mp_gMCS_I_gMCS_GMCS_LIST_TASKS",
                        #     checkboxGroupInput(inputId = "I_gMCS_GMCS_LIST_TASKS", label = ' ' ,
                        #                        choices = metTasks$EssentialTasks_CultureMedium$DESCRIPTION,
                        #                        selected = metTasks$EssentialTasks_CultureMedium$DESCRIPTION),
                        #     
                        #     br(),
                        # ),
                        
                        # shinyTree("I_gMCS_GMCS_LIST_TASKS_tree", checkbox = TRUE),
                        
                        
                        div(id = "mp_gMCS_I_gMCS_GMCS_LIST_TASKS_tree",
                            
                            br(),
                            shinyTree("I_gMCS_GMCS_LIST_TASKS_tree", checkbox = TRUE, themeIcons = FALSE,
                                      multiple = T,
                                      themeDots = TRUE, theme = 'proton', wholerow = F), # default, default-dark, proton
                            tags$head(tags$style(type = "text/css",".jstree-proton .jstree-clicked { background: #0000; color: #000; border-radius: 3px; box-shadow: inset 0 0 1px #0000; }")),
                            # tags$head(tags$style(type = "text/css",".jstree-proton .jstree-clicked { background: #0000; color: #000; border-radius: 3px; box-shadow: inset 0 0 1px; }")),
                            # tags$head(tags$style(type = "text/css",".jstree-proton .jstree-clicked { background: #0000; border-radius: 3px; box-shadow: inset 0 0 1px #0000; }")),
                            tags$head(tags$style(type = "text/css",".jstree-proton .jstree-anchor {margin: 1px 0 10px; height: auto !important;}")),
                            tags$head(tags$style(type = "text/css",".jstree-leaf {height: auto;}")),
                            tags$head(tags$style(type = "text/css",".jstree-leaf a {height: auto !important;}")),
                            tags$head(tags$style(type = "text/css",".jstree-proton a {white-space:normal !important; height: auto; }")),
                            tags$head(tags$style(type = "text/css",".jstree-proton li > ins {vertical-align:top; }")),
                        ),

                        div(id = "mp_gMCS_GMCS_LIST_max_length",
                            br(),
                            sliderInput(inputId = "I_gMCS_GMCS_LIST_max_length", "Maximum length of selected gMCSs:",
                                        min = 0, max = 40, value = c(40))
                        ),
                        
                        
                        # actionButton("I_5_plot_button",
                        #              "Plot!",
                        #              width = '100%',heigth = '100%',
                        #              style="color: #ffffff; background-color: #333333; border-color: #ffffff")
                        
           ),
           # Main panel
           mainPanel(# co_FDR = 0.4, co_log2FC = 2
             
             div(id = "mp_gMCS_3",
                 h3("Summary table of gMCS database and metabolic tasks:"),
                 br(),
                 div(withSpinner(DT::dataTableOutput(outputId = "O_gMCS_metabolic_tasks")), style = "font-size: 95%; width: 100%")
             ),
             div(id = "mp_gMCS_download",
                 fluidRow(
                   column(width = 5, downloadButton("O_gMCS_database_download","Download Rdata", width = '100%',heigth = '90%')),
                   column(width = 2),
                   column(width = 5, actionButton(inputId = "I_show_summary_images_gMCS_database","Show/hide summary images", width = '100%',heigth = '90%'))
                 ),
                 br(),
             ),
             
             div(id = "mp_gMCS_1",
                 h3("Summary images of gMCS database:"),
                 br(),
                 fluidRow(
                   column(width = 7, 
                          div(withSpinner(plotOutput(outputId = "O_gMCS_piechart_gmcs", height="600px", width="700px")
                          ), style = "font-size: 95%; width: 100%")
                   ),
                   column(width = 5, 
                          div(withSpinner(plotOutput(outputId = "O_gMCS_barchart_merged_gmcs", height="600px", width="500px")
                          ), style = "font-size: 95%; width: 100%")
                   )
                 )),
             div(id = "mp_gMCS_2",
                 br(),
                 div(withSpinner(plotOutput(outputId = "O_gMCS_barchart_gmcs", height="600px", width="1200px")
                 ), style = "font-size: 95%; width: 100%")
             ),
             
           )
         )
)


