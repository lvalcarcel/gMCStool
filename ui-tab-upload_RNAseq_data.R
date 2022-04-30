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


tabPanel(title = "2. Upload RNA-seq data",
         
         sidebarLayout(
           # Sidebar
           sidebarPanel(width = 3,
                        h3("Here you can upload the data"),
                        # h5(":"),
                        # h6(textOutput("O_txtout3")),
                        # h5("Number of samples:"),
                        # h6(textOutput("O_nCL_sel3")),
                        # h5("Number of Protein Targets:"),
                        # h6(textOutput("O_txtout_nEsGn3")),
                        # hr(style="border-color: black; background-color: black;color: black;"),
                        
                        
                        # Select list of gMCSs to use in the study
                        
                        # h4("Select the method to update the gene information:"),
                        radioButtons(inputId = "I_DataInput_mode", label = "Select the method to update the gene information:",
                                     choiceValues = list("text",
                                                         "tximport",
                                                         "Rdata"),
                                     choiceNames =  list("Text files",
                                                         "tximport files",
                                                         "Rdata from previous sessions"),
                                     selected = "text"),
                        
                        actionButton("I_Data_input_examples",
                                     "Show Examples",
                                     width = '100%',heigth = '100%',
                                     style="color: #ffffff; background-color: #333333; border-color: #ffffff")
           ),
           # Main panel
           mainPanel(# 
             h3("Data upload:"), 
             div(id = "mp2_1",
                 h4("File:"),
                 fileInput("Data_input_text_expr", "Choose text file which contains the gene expression (it will eliminate sample classification)",
                           width = "100%", accept = c(".txt",".tsv",".csv")
                 ), 
                 fileInput("Data_input_text_class", "Choose text file which contains sample information (it will overwrite changes in sample classification)",
                           width = "100%", accept = c(".txt",".tsv",".csv")
                 )),
             div(id = "mp2_2",
                 h4("tximport output:"),
                 fileInput("Data_input_tximport", "Choose .RData or .RDS file to upload tximport output (it will eliminate sample classification)",
                           width = "100%", accept = c(".Rdata", ".RDS")
                 ),
                 fileInput("Data_input_text_class_2", "Choose text file which contains sample information (it will overwrite changes in sample classification)",
                           width = "100%", accept = c(".txt",".tsv",".csv")
                 )),
             div(id = "mp2_3",
                 h4("Previous sessions and examples:"),
                 fileInput("Data_input_rdata", "Choose .RData file to upload previously saved results",
                           width = "100%", accept = c(".Rdata")
                 )),
             div(id = "mp2_4",
                 actionButton("I_Data_input_example_TPM",
                              "Load Example Data of Bcell and MM [TPM]",
                              width = '100%',heigth = '100%',
                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                 actionButton("I_Data_input_example_log2TPM",
                              "Load Example Data of Bcell and MM [log2(TPM+1)]",
                              width = '100%',heigth = '100%',
                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                 actionButton("I_Data_input_example_DepMap",
                              "Load Example Data of DepMap Cell Lines [TPM]",
                              width = '100%',heigth = '100%',
                              style="color: #ffffff; background-color: #333333; border-color: #ffffff"), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
             ),
             br(),
             div(id = "mp2",
                 h3("Summary of the input data:"),
                 h5("You can set here all the values for the sample class and cohort:"),
                 br(),
                 
                 fluidRow(column(width = 2, h5("Sample class summary:")), 
                          column(width = 10, div(withSpinner(tableOutput(outputId = "O_table_sample_class_cohort_1")),style = "font-size: 95%; width: 66%"))),
                 fluidRow(column(width = 2, h5("Sample cohort summary:")),
                          column(width = 10,  div(withSpinner(tableOutput(outputId = "O_table_sample_class_cohort_2")), style = "font-size: 95%; width: 66%"),)),
                 
                 
                 
                 
                 
                 fluidRow(column(width = 3, actionButton("I_Data_input_alpha_sort_class",
                                                         "Alphabetically sort sample class",
                                                         width = '100%',heigth = '90%',
                                                         style="color: #ffffff; background-color: #333333; border-color: #ffffff")),
                          column(width = 3, actionButton("I_Data_input_num_sort_class",
                                                         "Numerically sort sample cohort",
                                                         width = '100%',heigth = '90%',
                                                         style="color: #ffffff; background-color: #333333; border-color: #ffffff")),
                          # column(width = 2),
                          column(width = 3, actionButton("I_Data_input_alpha_sort_cohort",
                                                         "Alphabetically sort sample cohort",
                                                         width = '100%',heigth = '90%',
                                                         style="color: #ffffff; background-color: #333333; border-color: #ffffff")),
                          column(width = 3, actionButton("I_Data_input_num_sort_cohort",
                                                         "Numerically sort sample cohort",
                                                         width = '100%',heigth = '90%',
                                                         style="color: #ffffff; background-color: #333333; border-color: #ffffff"))
                 ),
                 br(),
                 div(withSpinner(DT::dataTableOutput(outputId = "O_table_sample_class_cohort")), style = "font-size: 95%; width: 100%"),
                 br(),
                 br(),
                 fluidRow(
                   column(width = 4,  downloadButton("O_data_input_download_1","Download gene expression (txt)", width = '100%')),
                   column(width = 4,  downloadButton("O_data_input_download_2","Download sample classification (txt)", width = '100%')),
                   column(width = 4,  downloadButton("O_data_input_download_rdata","Download gene expression and sample classification (RData", width = '100%'))
                 ),
                 # br(),
                 br(),
                 actionButton("I_mp2_next", "Next!", width = '100%',heigth = '200%',
                              style=STYLE_NEXT_BUTTONS), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                 
                 br(),
             ),
             div(id = "mp2_loading",
                 h3("Loading the data:"),
                 br(),
                 fluidRow(column(width = 2),
                          column(width = 6, img(src="loading-74.gif", align = "center",height='350px',width='350px')),
                          column(width = 2)),
                 br()
             )
             
           )
         )
)


