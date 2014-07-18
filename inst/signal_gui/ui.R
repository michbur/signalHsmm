library(shiny)

inputTextarea <- function(inputId, value="", nrows, ncols) {
  tagList(
    singleton(tags$head(tags$script(src = "textarea.js"))),
    tags$textarea(id = inputId,
                  class = "inputtextarea",
                  rows = nrows,
                  cols = ncols,
                  as.character(value))
  )
}


shinyUI(pageWithSidebar(

  headerPanel("Signal.hsmm"),
  
  sidebarPanel(
    fileInput('seq_file', 'Choose .fasta or .txt file',
              accept=c('text/plain'))
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Data input", verbatimTextOutput("summary")),
      tabPanel("Short output", tableOutput("pred_table")),
      tabPanel("Long output", verbatimTextOutput("pred_long"))
    )
  )
))