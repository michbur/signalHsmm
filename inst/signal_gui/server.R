library(shiny)
library(signal.hsmm)
library(shinyAce)
options(shiny.maxRequestSize=10*1024^2)



shinyServer(function(input, output) {
  
  prediction <- reactive({
    
    if (!is.null(input[["seq_file"]]))
      res <- read_txt(input[["seq_file"]][["datapath"]])
    input[["use_area"]]
    isolate({if (input[["text_area"]] != "")
      res <- read_txt(textConnection(input[["text_area"]]))
    })
    
    if(length(res) > 300) {
      #dummy error, just to stop further processing
      stop("Too many sequences.")
    } else {
      run_signal.hsmm(res)
    }
  })
  
  
  output$dynamic_ui <- renderUI({
    if(class(try(prediction(), silent = TRUE)) == "try-error") {
      div(tags$p("Waiting for valid input"),
          actionButton("use_area", "Submit data from field on right..."),
          fileInput('seq_file', '...or choose .fasta or .txt file:'),
          tags$p("Queries bigger than 300 sequences will be not processed. 
                 Use batch mode instead."))
    } else {
      div(tags$p("Be patient - bigger calculations take few minutes."),
          downloadButton("download_short", "Download short output"),
          downloadButton("download_long", "Download long output (without graphics)"),
          downloadButton("download_long_graph", "Download long output (with graphics)"),
          tags$p("Refresh the page to start a new query with signal.hsmm."))
    }
  })
  
  output$dynamic_panel <- renderUI({
    if(class(try(prediction(), silent = TRUE)) == "try-error") {
      aceEditor("text_area", value="", height = 150)
    } else {
      verbatimTextOutput("summary")
    }
  })
  
  
  output$pred_table <- renderTable({
    if(class(try(prediction(), silent = TRUE)) == "try-error") {
      data.frame(sp.probability = "No", 
                 sp.start = "sequence",
                 sp.end = "chosen")
    } else {
      pred2df(prediction())
    }
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
    if(class(try(prediction(), silent = TRUE)) == "try-error") {
      verbatimTextOutput("pred_long_null")
    } else {
      uiOutput("long_preds")
    }
  })
  
  output$pred_long_null <- renderPrint({
    cat("No sequence chosen.")
  })
  
  output$download_short <- downloadHandler(
    filename  = function() { 
      part_name <- strsplit(input[["seq_file"]][["name"]], ".", fixed = TRUE)[[1]][1]
      paste0(part_name, "_pred.csv") 
    },
    content <- function(file) {
      write.csv(pred2df(prediction()), file)
    }
  )
  
  output$download_long <- downloadHandler(
    filename  = function() { 
      part_name <- strsplit(input[["seq_file"]][["name"]], ".", fixed = TRUE)[[1]][1]
      paste0(part_name, "_pred.txt") 
    },
    content <- function(file) {
      sink(file, type = "output")
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
      part_name <- strsplit(input[["seq_file"]][["name"]], ".", fixed = TRUE)[[1]][1]
      paste0(part_name, "_pred.html") 
    },
    content <- function(file) {
      knitr:::knit(input = "signalhsmm_report.Rmd", 
                   output = "signalhsmm_report.md", quiet = TRUE)
      markdown:::markdownToHTML("signalhsmm_report.md", file)
    }
  )
})
