require(shiny)
shinyUI(pageWithSidebar(
  headerPanel("Phytobenthos Classification"),
  sidebarPanel(
    fileInput('file1', 'Choose CSV File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    tags$hr(),
    checkboxInput('lake', 'Lake', TRUE),
    checkboxInput('river', 'River', TRUE)            
  ),
  mainPanel(
    tableOutput('contents')
  )
))
