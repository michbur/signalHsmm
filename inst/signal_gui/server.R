library(shiny)
library(signal.hsmm)
options(shiny.maxRequestSize=15*1024^2)

shinyServer(function(input, output) {
  
  prediction <- reactive({
    if(grepl("fasta", input[["seq_file"]][["name"]])) {
      input_dat <- seqinr:::read.fasta(input[["seq_file"]][["datapath"]])
    } 
    #todo - read from .txt files
    run_signal.hsmm(input_dat)
  })
  
  
  
  output$pred_table <- renderTable({
    if(is.null(input[["seq_file"]])) {
      data.frame(sp.probability = "No", 
                 sp.start = "file",
                 sp.end = "chosen")
    } else {
      pred2df(prediction())
    }
  })
  
  output$summary <- renderPrint({
    if(is.null(input[["seq_file"]])) {
      cat("No file chosen.")
    } else {
      summary(prediction())
    }
  })
  
  output$pred_long <- renderPrint({
    if(is.null(input[["seq_file"]])) {
      cat("No file chosen.")
    } else {
      for (i in 1L:length(prediction())) {
        cat(names(prediction())[i])
        cat("\n\n")
        summary(prediction()[[i]])
        cat("\n\n")
      }
    }
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
