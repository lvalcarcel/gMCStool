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



tabPanel("Help",
         fluidPage(
           includeMarkdown("instructions/gettingStarted_help_v4.Rmd"),
           # includeHTML("instructions/gettingStarted_help_v3.knit.html"),
           br(),
           br(),
           fluidRow(column(width = 3),
                    column(width = 6,
                           downloadButton("O_help_download","Download section in PDF", width = '100%',heigth = '90%'))
                    
           )
         ))

# div(img(src = "welcome1.png", width = 900), style="text-align: center;")),