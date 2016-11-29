#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(markdown)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    
    # Application title
    titlePanel("MDCK RNA-Seq analysis - Version 0.1"),

    mainPanel(
      tabsetPanel(type = "tabs",
        tabPanel("Volcano plot - EMT markers - TGFb vs WT", plotOutput("volcano_tgfb")),
        tabPanel("EMT Markers - TGFb vs WT", DT::dataTableOutput("emtMarker_de_tgfb")),
        tabPanel("Volcano plot - EMT markers - shZ KD vs WT", plotOutput("volcano_shZ")),
        tabPanel("EMT Markers - shZ KD vs WT", DT::dataTableOutput("emtMarker_de_shZ")),
        tabPanel("Global DE - TGFb vs WT", DT::dataTableOutput("global_de_tgfb")),
        tabPanel("Global DE - shZ KD vs WT", DT::dataTableOutput("global_de_shZ"))
    )
  )
))