library(shiny)



shinyUI(fluidPage(

  headerPanel("Signal.hsmm"),
  
  sidebarPanel(
    tags$div(class="header", checked=NA,
             tags$p("signal.hsmm detects signal peptides in eukaryotic proteins.")),
    tags$div(class="header", checked=NA,
             tags$p(" ")),
    fileInput('seq_file', 'Choose .fasta or .txt file:',
              accept=c('text/plain')),
    downloadButton("download_short", "Download short output"),
    downloadButton("download_long", "Download long output")

  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Data input", verbatimTextOutput("summary")),
      tabPanel("Short output", tableOutput("pred_table")),
      tabPanel("Long output", verbatimTextOutput("pred_long"))
    )
  )
))
