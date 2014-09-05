library(shiny)
#library(shinyAce)



shinyUI(fluidPage(
  
  headerPanel("signal.hsmm"),
  
  sidebarLayout(
    sidebarPanel(
      uiOutput("dynamic_ui"),
      includeMarkdown("readme.md")
    ),
    
    mainPanel(
      uiOutput("dynamic_tabset")    
    )
  )))
