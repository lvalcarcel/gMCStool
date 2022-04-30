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



tabPanel("Overview",
         useShinyjs(),  # Set up shinyjs
         
         fluidPage(
           includeMarkdown("instructions/gettingStarted_overview_v2.Rmd"),
           
           fluidRow(column(width = 2),
                    column(width = 8,
                           actionButton("I_mp0_next", "Start!", width = '100%',heigth = '200%',
                                        style=STYLE_NEXT_BUTTONS), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                    )
           ),
           br(),
           fluidRow(column(width = 2),
                    column(width = 8,
                           actionButton("I_mp0_next_to_download", "Go to gMCS database", width = '100%',heigth = '100%',
                                        style=STYLE_NEXT_BUTTONS), # Look for text redistribution inside the actionButton !!!!!!!!!!!!!
                           
                    )
           )
         )
         
         
)



# div(img(src = "welcome1.png", width = 900), style="text-align: center;")),