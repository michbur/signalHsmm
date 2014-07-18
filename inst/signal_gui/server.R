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
    pred2df(prediction())
  })
  
  output$summary <- renderPrint({
    if(is.null(input[["seq_file"]])) {
      cat("No file chosen.")
    } else {
      summary(prediction())
    }
  })
  
  output$pred_long <- renderPrint({
    for (i in 1L:length(prediction())){
      cat(names(prediction())[i])
      cat("\n\n")
      summary(prediction()[[i]])
      cat("\n\n")
    }
  })
  
})
