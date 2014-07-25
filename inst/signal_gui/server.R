library(shiny)
library(signal.hsmm)
library(shinyAce)
options(shiny.maxRequestSize=10*1024^2)

shinyServer(function(input, output) {
  
  
  
  input_file <- reactive({
    res <- read_txt(input[["seq_file"]][["datapath"]])
    #     input[["use_area"]]
    #     isolate(input[["text_area"]])
    res
  })
  
  prediction <- reactive({
    #     print(input_area())
    #     print()
    #     inputs <- list(input_area(), input_file())
    
    run_signal.hsmm(input_file())
  })
  
  
  output$dynamic_ui <- renderUI({
    if (is.null(input[["seq_file"]])) {
      div(actionButton("use_area", "Submit data from field on right..."),
          fileInput('seq_file', '...or choose .fasta or .txt file:'))
    } else {
      div(downloadButton("download_short", "Download short output"),
          downloadButton("download_long", "Download long output (no graphics)"))
    }
  })
  
  output$dynamic_panel <- renderUI({
    if (is.null(input[["seq_file"]])) {
      aceEditor("text_area", value="", height = 150)
    } else {
      verbatimTextOutput("summary")
    }
  })
  
  
  output$pred_table <- renderTable({
    if(is.null(input[["seq_file"]])) {
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
  
  
  #   output$pred_long <- renderPrint({
  #     if(is.null(input[["seq_file"]])) {
  #       cat("No sequence chosen.")
  #     } else {
  #       for (i in 1L:length(prediction())) {
  #         cat(names(prediction())[i])
  #         cat("\n\n")
  #         summary(prediction()[[i]])
  #         cat("\n\n")
  #       }
  #     }
  #   })
  
  output$long_preds <- renderUI({
    long_preds_list <- lapply(1L:length(prediction()), function(i) {
      list(plotOutput(paste0("plot", i)), verbatimTextOutput(paste0("summ", i)))
    })
    do.call(tagList, unlist(long_preds_list, recursive = FALSE))
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  for (i in 1L:5) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      
      output[[paste0("plot", my_i)]] <- renderPlot(plot(prediction()[[my_i]]))
      output[[paste0("summ", my_i)]] <- renderPrint(summary(prediction()[[my_i]]))
    })
  }
  
  
  output$pred_long <- renderUI({
    if(is.null(input[["seq_file"]])) {
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
        cat(names(prediction())[i])
        cat("\n\n")
        summary(prediction()[[i]])
        cat("\n\n")
      }
      sink()
    }
  )
  
})
