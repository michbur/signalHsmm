library(shiny)
#library(shinyAce)



shinyUI(fluidPage(
  
  headerPanel("Signal.hsmm"),
  
  sidebarLayout(
  sidebarPanel(
    uiOutput("dynamic_ui"),
    includeMarkdown("readme.md")
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Data input", uiOutput("dynamic_panel")),
      tabPanel("Short output", tableOutput("pred_table")),
      tabPanel("Long output (with graphics)", uiOutput("pred_long"))
    )
  )
)))
