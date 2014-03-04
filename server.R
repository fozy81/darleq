require(shiny)
shinyServer(function(input, output) {
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects and uploads a 
    # file, it will be a data frame with 'name', 'size', 'type', and 'datapath' 
    # columns. The 'datapath' column will contain the local filenames where the 
    # data can be found.
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    testdata <-read.csv(inFile$datapath)
    source("darleqFunc.R")
    dataTDI <- darleqFunc(testdata)
    dataTDI <- dataTDI[,order(names(dataTDI), decreasing = TRUE)]
    row.names(dataTDI) <- NULL  # remove row names not required for display
    dataTDI$'SAMPLE ID' <- as.character(floor(as.numeric(dataTDI$'SAMPLE ID'))) # round sampleID
    lake <- input$lake
    if (input$lake == TRUE & input$river == FALSE) # return different bits of table depending on river or lake
      return(dataTDI[,grepl("LAKE*|SAMPLE*", names(dataTDI))])
    if (input$river == TRUE & input$lake == FALSE)
      return(dataTDI[,grepl("RIVER*|SAMPLE*", names(dataTDI))])
    
    #stats[grepl("CD.Mean*", rownames(stats)), ]
    #df[,grep('^x',names(df))] 
    #df[ , -which(names(df) %in% c("z","u"))]
    dataTDI
  })
})
