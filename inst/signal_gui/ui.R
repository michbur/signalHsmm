library(shiny)

shinyUI(fluidPage(
  
  headerPanel("signalHsmm"),
  
  sidebarLayout(
    sidebarPanel(
      includeMarkdown("readme.md"),
      pre(includeText("prots.txt")),
      uiOutput("dynamic_ui")
    ),
    
    mainPanel(
      uiOutput("dynamic_tabset")    
    )
  )))
