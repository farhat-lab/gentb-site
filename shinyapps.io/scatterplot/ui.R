library(shiny)

library(jsonlite)
library(tidyr)
library(rCharts)

##----- Drug list

drug_list <- c('inh','rif','pza', 'emb','str','eth','kan', 'cap',
               'amk', 'cip', 'levo', 'oflx', 'pas') 

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Scatterplot"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput('selected_strain', 'Choose Strain', c("01-R1466", "01-R1467", "01-R1469")),
      selectInput('selected_drug', 'Choose Drug', drug_list),
      downloadButton('cars', 'Download data')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))