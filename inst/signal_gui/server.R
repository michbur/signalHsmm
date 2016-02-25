library(shiny)
library(signalHsmm)
library(shinyAce)
options(shiny.maxRequestSize=10*1024^2)



shinyServer(function(input, output) {
  
  prediction <- reactive({
    
    if (!is.null(input[["seq_file"]]))
      res <- read_txt(input[["seq_file"]][["datapath"]])
    input[["use_area"]]
    isolate({
      if (!is.null(input[["text_area"]]))
        if(input[["text_area"]] != "")
          res <- read_txt(textConnection(input[["text_area"]]))
    })
    
    if(exists("res")) {
      if(length(res) > 300) {
        #dummy error, just to stop further processing
        stop("Too many sequences.")
      } else {
        run_signalHsmm(res)
      }
    } else {
      NULL
    }
  })
  
  
  output$dynamic_ui <- renderUI({
    if(!is.null(prediction())) {
      div(tags$h3("Download results"),
          tags$p(""),
          downloadButton("download_short", "Download short output"),
          downloadButton("download_long", "Download long output (without graphics)"),
          downloadButton("download_long_graph", "Download long output (with graphics)"),
          tags$p("Refresh page (press F5) to start a new query with signalHsmm."))
    }
  })
  
  
  output$pred_table <- renderTable({
    pred2df(prediction())
  })
  
  output$summary <- renderPrint({
    summary(prediction())
  })
  
  
  output$long_preds <- renderUI({
    long_preds_list <- lapply(1L:length(prediction()), function(i) {
      list(plotOutput(paste0("plot", i)), verbatimTextOutput(paste0("summ", i)))
    })
    do.call(tagList, unlist(long_preds_list, recursive = FALSE))
  })
  
  
  for (i in 1L:300) {
    local({
      my_i <- i
      
      output[[paste0("plot", my_i)]] <- renderPlot(plot(prediction()[[my_i]]))
      output[[paste0("summ", my_i)]] <- renderPrint(summary(prediction()[[my_i]]))
    })
  }
  
  
  output$pred_long <- renderUI({
    uiOutput("long_preds")
  })
  
  output$dynamic_tabset <- renderUI({
    if(is.null(prediction())) {
      
      tabPanel(title = "Sequence input",
               h3("Paste sequences (FASTA format required) into the field below:"), 
               tags$style(type="text/css", "textarea {width:100%}"),
               tags$textarea(id = "text_area", rows = 22, cols = 60, ""),
               p(""),
               actionButton("use_area", "Submit data from field above"),
               p(""),
               fileInput('seq_file', 'Submit .fasta or .txt file:'))
      
      
    } else {
      tabsetPanel(
        tabPanel("Input summary", verbatimTextOutput("summary")),
        tabPanel("Short output", tableOutput("pred_table")),
        tabPanel("Long output (with graphics)", uiOutput("pred_long"))
      )
    }
  })
  
  #name for downloads
  file_name <- reactive({
    if(is.null(input[["seq_file"]][["name"]])) {
      part_name <- "signalhsmm_results"
    } else {
      part_name <- strsplit(input[["seq_file"]][["name"]], ".", fixed = TRUE)[[1]][1]
    }
    part_name
  })
  
  
  output$download_short <- downloadHandler(
    filename  = function() { 
      paste0(file_name(), "_pred.csv") 
    },
    content <- function(file) {
      write.csv(pred2df(prediction()), file)
    }
  )
  
  output$download_long <- downloadHandler(
    filename  = function() { 
      paste0(file_name(), "_pred.txt") 
    },
    content <- function(file) {
      sink(file, type = "output")
      cat("Input file name: ", ifelse(is.null(input[["seq_file"]][["name"]]), "none",
                                      input[["seq_file"]][["name"]]), "\n\n")
      cat(paste0("Date: ", Sys.time()), "\n\n")
      for (i in 1L:length(prediction())) {
        cat("\n\n")
        summary(prediction()[[i]])
        cat("\n\n")
      }
      sink()
    }
  )
  
  output$download_long_graph <- downloadHandler(
    filename  = function() { 
      paste0(file_name(), "_pred.html") 
    },
    content <- function(file) {
      knitr:::knit(input = "signalhsmm_report.Rmd", 
                   output = "signalhsmm_report.md", quiet = TRUE)
      markdown:::markdownToHTML("signalhsmm_report.md", file)
    }
  )
})