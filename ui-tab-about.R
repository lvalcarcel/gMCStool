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


tabPanel("About",
         
         
         flipBox(
           id = "myflipbox",
           trigger = "hover",
           # main_img = "tecnunLogoBigIcon.jpg",
           # header_img = 'ibaeta.jpg',
           # front_title = "Biomedical Engineering and Sciences Department. University of Navarra - Tecnun.",
           # back_title = "Contact Info",
           width = 12,
           front = div(
             class = "text-center",
             h1("Biomedical Engineering and Sciences Department. Computational Biology group at Tecnun- University of Navarra."),
            
             
             HTML(paste("<br>Our research group focus on development and application of computational systems biology approaches to analyse -omics data.
                        This tool was developed in close collaboration of Onco-hematology group at ",
                        a('CIMA- University of Navarra', 
                          href='https://cima.cun.es/investigacion/programas-investigacion/programa-investigacion-hemato-oncologia')),
                  ".<br>"),
             
             
             br(),     
             hr(style="border-color: #932636; background-color: #932636;color: #932636; width: 50%"),
             h3("gMCStool Team Members"),
             br(),
             fluidRow(
               userBox(
                 title = userDescription(
                   title = "Luis V. Valc\U00E1rcel",
                   subtitle = HTML("<b>PhD Student</b><br>Biomedical Engineering"),
                   type = NULL,
                   image = "valcarcel_luis.jpg"
                 ),
                 width = 4,
                 color = "aqua-active",
                 closable = FALSE,
                 collapsible = FALSE,
                 footer = ""
               ),
               userBox(
                 title = userDescription(
                   title = "I\U00F1igo Apaolaza",
                   subtitle = HTML("<b>Postdoc</b><br>Applied Engineering PhD"),
                   type = NULL,
                   image = "apaolaza_emparanza_inigo.jpg"
                 ),
                 width = 4,
                 color = "aqua-active",
                 closable = FALSE,
                 collapsible = FALSE,
                 footer = ""
               ),
               userBox(
                 title = userDescription(
                   title = "Francisco Javier Planes Pedre\U00F1o",
                   subtitle = HTML("<b>TECNUN Deputy Director of Research</b><br>Industrial Engineer PhD"),
                   type = NULL,
                   image = "planes_pedreno_francis.jpg"
                 ),
                 width = 4,
                 color = "aqua-active",
                 closable = FALSE,
                 collapsible = FALSE,
                 footer = ""
               )
             ),
             
             br(),br(),br(),
             fluidRow(
               userBox(
                 title = userDescription(
                   title = "Edurne San Jose",
                   subtitle = HTML("<b>Researcher</b><br>Biology PhD"),
                   type = NULL,
                   image = "edurne-san-jose-eneriz.jpg"
                 ),
                 width = 4,
                 color = "aqua-active",
                 closable = FALSE,
                 collapsible = FALSE,
                 footer = ""
               ),
               userBox(
                 title = userDescription(
                   title = "Xabier Agirre Ena",
                   subtitle = HTML("<b>Researcher</b><br>Biology PhD"),
                   type = NULL,
                   image = "agirre_xabier.jpg"
                 ),
                 width = 4,
                 color = "aqua-active",
                 closable = FALSE,
                 collapsible = FALSE,
                 footer = ""
               ),
               userBox(
                 title = userDescription(
                   title = "Felipe Pr\U00F3sper Cardoso",
                   subtitle = HTML("<b>CIMA Director of Cellular Therapy</b><br>Medical Doctor"),
                   type = NULL,
                   image = "dr-felipe-prosper_2.jpg"
                 ),
                 width = 4,
                 color = "aqua-active",
                 closable = FALSE,
                 collapsible = FALSE,
                 footer = ""
               )
               
             ),
             br()
             
             # actionButton("toggle", "Contact")
             
           ), # end front
           
           # front_btn_text = "Contact",
           # back_btn_text = "Back",
           # back_content = tagList(
           back = div(
             
             br(),br(),
             class = "text-center",
             h1("Contact Info"),
             
             fluidRow(
               shinydashboard::box(
                 solidHeader = FALSE,
                 title = "TECNUN - Miram\U00F3n Campus",
                 background = NULL,
                 width = 12,
                 status = "success",
                 footer = fluidRow(
                   br(),
                   column(
                     width = 6,
                     HTML(paste("Email: fplanes@tecnun.es"))
                   ),
                   column(
                     width = 6,
                     HTML(paste("Tel. 943 21 98 77"))
                   )
                 )
               ),
             ),
             br(),br(),
             fluidRow(
               htmlOutput("iframe")
             )#,
             # actionButton("toggle", "Info")
             
             
           ) # end back
           
         )
         
         
         
)

# div(img(src = "welcome1.png", width = 900), style="text-align: center;")),