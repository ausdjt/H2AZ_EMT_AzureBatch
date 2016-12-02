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
library(DT)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    
    # Application title
    titlePanel("MDCK RNA-Seq analysis - Version 0.2"),

    mainPanel(
      tabsetPanel(type = "tabs",
        tabPanel("Volcano plot - EMT markers - TGFb vs WT", plotOutput("volcano_tgfb")),
        tabPanel("EMT Markers - TGFb vs WT", DT::dataTableOutput("emtMarker_de_tgfb")),
        tabPanel("Volcano plot - EMT markers - shZ KD vs WT", plotOutput("volcano_shZ")),
        tabPanel("EMT Markers - shZ KD vs WT", DT::dataTableOutput("emtMarker_de_shZ")),
        tabPanel("Global DE - TGFb vs WT", DT::dataTableOutput("global_de_tgfb")),
        tabPanel("Global DE - shZ KD vs WT", DT::dataTableOutput("global_de_shZ")),
        tabPanel("Quality Control report", htmlOutput("inc_html")),
        tabPanel("Interactive volcano plot - EMT signature genes", 
                 fluidRow(
                    column(width=3,
                           div(class = "option-group",
                               radioButtons("dataset", "Data set",
                                            choices = c("conditionMDCKshZ", "conditionMDCKTGFb"), inline = FALSE),
                               selectInput("highlight", "Gene group to highlight",
                                          c("none" = "none",
                                            "Epithelial markers" = "epi",
                                            "Mesenchymal markers" = "mes",
                                            "Both" = "both"),
                                          selected = "none",
                                          selectize = FALSE)
                              )),
                    column(width = 9, uiOutput("plotui"))
                    ), #fluidRow 1
                 fluidRow(DT::dataTableOutput("plot_brushed_points"))
        )
      )
    )
  )
)