library(leaflet)
library(shiny)
library(ggvis)
library(ggplot2)

shinyUI(
  bootstrapPage(

    tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css')),

    titlePanel("Select a country on the map"),

    sidebarLayout(

      mainPanel(

        tabsetPanel(

          tabPanel("Drug resistance",
                   h4(textOutput("country2")),
                   downloadButton('download_country_data', 'Download country data'),
                   downloadButton('download_world_data', 'Download world data'),
                   h4(''),
                   ggvisOutput("resistance_plot"),
                   uiOutput("resistance_plot_ui")
          ),

          tabPanel("Lineage distribution",
                   h4(textOutput("country1")),
                   h4(''),
                   textOutput("lineage_indicator"),
                   conditionalPanel(condition="output.lineage_indicator == ' '", ggvisOutput("lineage_plot")),
                   conditionalPanel(condition="output.lineage_indicator == ' '", uiOutput("lineage_plot_ui")),
                   conditionalPanel(condition="output.lineage_indicator == '  '", textOutput("error_message"))
          ),

          tabPanel("Mutations",
                   h4(textOutput("country3")),
                   dataTableOutput('mutation_table')
          ),

          tabPanel("Map setup",
                   checkboxInput("resize", "Reduce circles", FALSE),
                   checkboxInput("greyscale", "Grey scale", FALSE))
        )
      ),

      sidebarPanel(
        leafletOutput("map")
      ),
      position = "right"
    )
  )
)

