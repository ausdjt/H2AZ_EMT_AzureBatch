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
fluidPage(
  
  # Application title
  titlePanel("MDCK RNA-Seq analysis - Version 0.1"),
  
  fluidPage(
    includeHTML("mdck_rna-seq.html")
  )
)
