require(shiny)
shinyUI(pageWithSidebar(
  headerPanel("Phytobenthos Classification"),
  sidebarPanel(
    fileInput('file1', 'Choose CSV File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    tags$hr(),
    checkboxInput('lake', 'Lake', TRUE),
    checkboxInput('river', 'River', TRUE),  
    downloadButton('downloadTest','Dowload Output When Analysis Complete'),
   helpText(a("Download example demo test data here", href="https://raw.github.com/fozy81/darleq/master/testdata.csv", target="_blank")),
   helpText(a("Report issues or view the code for this site on Github", href="https://github.com/fozy81/darleq/issues", target="_blank"))
  ),
  mainPanel(
    tableOutput('table')
  )
))
