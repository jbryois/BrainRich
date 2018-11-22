library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Brain Cell Type Enrichment"),
  
  # Sidebars
  sidebarLayout(
    sidebarPanel(
      radioButtons('sep', 'Separator',
                   c(Tab='\t',
                     Space=' ',
                     Comma=',',
                     xls='xls',
                     xlsx='xlsx'),
                   ' '),
      checkboxInput('header', 'Header', TRUE),
      fileInput('file1', 'Input File'),
      uiOutput("choose_column"),
      uiOutput("top"),
      tableOutput("Number_genes"), 
      plotOutput("hist",height = "200px",width="75%"),
      actionButton("ewce_run", "Run EWCE"),
      numericInput("perm", "Permutations:", 99999),
      br(),
      downloadButton('save_plot',label = "Save Plot BrainRich"),
      downloadButton('save_table',label = "Save Table BrainRich"),
      br(),
      downloadButton('save_plot2',label = "Save Density Plot BrainRich"),
      br(),
      downloadButton('save_plot_ewce',label = "Save Plot EWCE"),
      downloadButton('save_table_ewce',label = "Save Table EWCE")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("Enrichment_p"),
      br(),
      #plotOutput("Cells"),
      #br(),
      plotOutput("EWCE_plot")
      #dataTableOutput("regression_table")
    )
  )
))


