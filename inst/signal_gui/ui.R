library(shiny)
#library(shinyAce)



shinyUI(fluidPage(
  
  headerPanel("Signal.hsmm"),
  
  sidebarLayout(
  sidebarPanel(
    tags$div(class="header", checked=NA,
             tags$p("signal.hsmm detects signal peptides in eukaryotic proteins.")),
    tags$p(" "),
    uiOutput("dynamic_ui")
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Data input", uiOutput("dynamic_panel")),
      tabPanel("Short output", tableOutput("pred_table")),
      tabPanel("Long output (with graphics)", uiOutput("pred_long"))
    )
  )
)))
