require(shiny)
shinyServer(function(input, output) {
 
  datasetInput <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    testdata <-read.csv(inFile$datapath)
    source("darleqFunc.R")
    dataTDI <- darleqFunc(testdata)
    dataTDI <- dataTDI[,order(names(dataTDI), decreasing = TRUE)]
    dataTDI$'SAMPLE ID' <- as.character(floor(as.numeric(dataTDI$'SAMPLE ID'))) # round sampleID
    row.names(dataTDI) <- NULL  # remove row names not required for display
    lake <- input$lake
    if (input$lake == TRUE & input$river == FALSE) # return different bits of table depending on river or lake
      return(dataTDI[,grepl("LAKE*|SAMPLE*", names(dataTDI))])
    if (input$river == TRUE & input$lake == FALSE)
      return(dataTDI[,grepl("RIVER*|SAMPLE*", names(dataTDI))])
    
    dataTDI
  })
  
  output$table <- renderTable({
    datasetInput()
  })
  
 
  output$downloadTest <- downloadHandler(
    filename = function() { paste(input$file1, '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file)
    }
  )

})
