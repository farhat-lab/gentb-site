library(shiny)
library(leaflet)
library(RColorBrewer)
library(maps)
library(data.table)
library(dplyr)
library(ggvis)
library(tidyr)
library(ggplot2)

# https://github.com/rstudio/shiny-examples/blob/master/063-superzip-example/server.R

gen_phen_country2 <- fread("dv_genotype_phenotype_with_country.csv")[, V1 := NULL]

gen_phen_country <- gen_phen_country2[, list(strain = strain,
              snpname = snpname,
              sub.ins.del = sub.ins.del,
              coding.silent = coding.silent,
              genomic.coordinate = genomic.coordinate,
              genetic.change = genetic.change,     
              locus = region,
              amino.acid.change = amino.acid.change,
              variant = coding.silent,
              # unclear,           
              code = Country_code,                  
              name = Name,               
              drtype = ResistanceType,             
              rinh = IsoniazidDST,                
              reth =EthionamideDST,  
              rrif = RifampicinDST,                 
              rrfb = RifabutinDST,
              remb = EthambutolDST,                 
              rpza = PyrazinamideDST,                 
              rstr = StreptomycinDST,
              ramk = AmikacinDST,        
              rcap = CapreomycinDST,                 
              rkan = KanamycinDST, 
              roflx = OfloxacinDST,
              rcip = CiprofloxacinDST,             
              rlevo = LevofloxacinDST,
              rmoxi = MoxifloxacinDST,
              rgati = GatifloxacinDST,                
              rpas = ParaAminoSalicylicAcidDST,
              rcys = CycloserineDST,     
              rtha = ThioacetazoneDST,                 
              rpro = ProthionamideDST,   
              rclof = ClofazimineDST,
              rclar = ClarithromycinDST,                
              ramoxclav = AmoxicillinClavulanateDST,
              rlin = LinezolidDST,
              date = Date,                 
              country = Country,        
              city = City,
              setting = Setting,              
              source = Source,
              ptage = PatientAge,
              ptsex = PatientSex,                
              hivstatus = HIVStatus,
              spfamily_parentstrain = Spfamily_parentstrain,
              RFLPfamily = RFLPfamily,           
              patientid = PatientID,
              spoligo_octal = Spoligo_octal,
              lat = lat,                  
              lon = lon)]

gen_phen_country <- data.frame(gen_phen_country)

table_gen_phen_country <- tbl_dt(gen_phen_country)

tidy_gen_phen_country <- gen_phen_country %>%
  gather(drug, resistance, rinh:rlin)

summary_gen_phen_country <- table_gen_phen_country %>%
  group_by(code) %>%
  summarize(
    country = country,
    resistant = sum(drtype %in% c("MDR", "XDR")),
    sub.ins.del = sub.ins.del,
    coding.silent = coding.silent,
    genomic.coordinate = genomic.coordinate,
    genetic.change = genetic.change,
    locus = locus,
    amino.acid.change = amino.acid.change,
    variant = variant
  )


gen_phen_country_by_strain2 <- fread("dv_genotype_phenotype_with_country_by_strain.csv", 
                                     na.strings = c("", " "))[, V1 := NULL]

gen_phen_country_by_strain <- gen_phen_country_by_strain2[, list(
              strain = strain,
              drtype = ResistanceType,             
              rinh = IsoniazidDST,                
              reth =EthionamideDST,  
              rrif = RifampicinDST,                 
              rrfb = RifabutinDST,
              remb = EthambutolDST,                 
              rpza = PyrazinamideDST,                 
              rstr = StreptomycinDST,
              ramk = AmikacinDST,        
              rcap = CapreomycinDST,                 
              rkan = KanamycinDST, 
              roflx = OfloxacinDST,
              rcip = CiprofloxacinDST,             
              rlevo = LevofloxacinDST,
              rmoxi = MoxifloxacinDST,
              rgati = GatifloxacinDST,                
              rpas = ParaAminoSalicylicAcidDST,
              rcys = CycloserineDST,     
              rtha = ThioacetazoneDST,                 
              rpro = ProthionamideDST,   
              rclof = ClofazimineDST,
              rclar = ClarithromycinDST,                
              ramoxclav = AmoxicillinClavulanateDST,
              rlin = LinezolidDST,
              country = Country,              
              code = Country_code,
              lat = lat,
              lon = lon,            
              spfamily_parentstrain = Spfamily_parentstrain,
              RFLPfamily = RFLPfamily)]

country_strain_data <- gen_phen_country_by_strain %>%
  group_by(code) %>%
  summarize(
    country = country,
    strain = n(),
    resistant = sum(drtype %in% c("MDR", "XDR")),
    lat = lat,
    lon = lon) %>% 
  distinct() %>%
  mutate(percent = 100 * resistant / strain)

country_strain_data <- as.data.frame(country_strain_data)

row.names(country_strain_data) <- country_strain_data$code

gen_phen_country_by_strain <- gen_phen_country_by_strain %>%
  mutate(family = ifelse(spfamily_parentstrain == "",
                         RFLPfamily, 
                         spfamily_parentstrain)) %>% data.frame()

country_lineage <- as.data.frame.matrix(table(gen_phen_country_by_strain$code, 
                                              gen_phen_country_by_strain$family))

shinyServer(function(input, output, session) {
  
  makeReactiveBinding("selected_country")
  
  marker_size <- reactive({
    ifelse(input$resize, 0.5, 1)
  })
  
  color_total_strains <- reactive({
    ifelse(input$greyscale, "grey", "#4A9")
  })
  
  color_resistant_strains <- reactive({
    ifelse(input$greyscale, "#403e3e", "#bb0707")
  })

  color_sensitive_strains <- reactive({
    ifelse(input$greyscale, "grey", "#1b443d")
  })
  
  color_intermediate_strains <- reactive({
    ifelse(input$greyscale, "white", "#b4ddd6")
  })

  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("Stamen.TonerLite") %>%
      setView(lng = 13, lat = 48, zoom = 1)
  })
  
  observe({
    if (is.null(input$map_click)){
      return()
    }
    selected_country <- NULL
  })
  
  observe({
    if (nrow(country_strain_data) == 0){
      return()
    }
    
    # Rescale values to be in [a, b]
    b <- 30
    a <- 5
    radii1 <- ifelse(country_strain_data$strain != 0,
                     (b - a) * 
                       (country_strain_data$strain  - min(country_strain_data$strain))  / 
                       (max(country_strain_data$strain) - min(country_strain_data$strain)) + a,
                     0) * input$map_zoom * marker_size()
    radii2 <- ifelse(country_strain_data$strain != 0, 
                     radii1 * 
                       country_strain_data$resistant / country_strain_data$strain, 
                     0)
    
    # Create unique layerIDs
    ids1 <- paste(country_strain_data$code, "strain", sep = ".")
    ids2 <- paste(country_strain_data$code, "resistant", sep = ".")
    
    # Popup Text
    popup_content <- paste(
      "<h4 style='border-bottom: thin dotted #43464C;
      padding-bottom:4px; margin-bottom:4px;
      font-family: Tahoma, Geneva, sans-serif;
      color:#43464C;'>", 
      country_strain_data$country, "</h4>
      <span style='color:#43464C;'>",
      "Number of total strains: ","</h4> <strong>",
      prettyNum(country_strain_data$strain, big.mark = ','), "<br/>",
      "</strong> <span style='color:#43464C;'>",
      "Number of resistant strains: ", "<strong>",
      prettyNum(country_strain_data$resistant, big.mark = ','), "</strong>")
    
    # Add circle markers and legend for Total and Resistant Strains
    leafletProxy("map") %>%
      clearShapes() %>%
      addCircleMarkers(lat = country_strain_data$lat,
                       lng = country_strain_data$lon,
                       layerId = ids1,
                       radius = radii1,
                       fill = TRUE,
                       fillColor = color_total_strains(),
                       color = color_total_strains(),
                       fillOpacity = .7,
                       opacity = .7,
                       weight = 1,
                       stroke = FALSE,
                       list(
                         clickable = TRUE,
                         riseOnHover = TRUE
                       ), popup = popup_content) %>%
      addCircleMarkers(lat = country_strain_data$lat,
                       lng = country_strain_data$lon,
                       layerId = ids2,
                       radius = radii2,
                       fill = TRUE,
                       fillColor = color_resistant_strains(),
                       color = color_resistant_strains(),
                       fillOpacity = 0.5,
                       opacity = 0.5,
                       weight = 1,
                       stroke = FALSE,
                       list(
                         clickable = TRUE,
                         riseOnHover = TRUE
                        ), popup = popup_content) %>%
      addLegend(
            position = 'bottomright',
            layerId = ids2,
            colors = c("#7CC4B9", "#983D35"),
            labels = c("Total Strains", "Resistant Strains"), opacity = 0.9)
    
  })
  
  observe({
    
    event <- input$map_marker_click
    
    if (is.null(event)){
      return()
    }

    isolate({
      country <- country_strain_data[paste(country_strain_data$code, 
                                           "strain", sep = ".") == event$id
                                    | paste(country_strain_data$code, 
                                           "resistant", sep = ".") == event$id, ]
      selected_country <- country

    })
    
# @asiegmann: we remove empty resistance types ... should we look for NA or
# actually include these? will this get all the bars to total to Total #? Play around with this, group by country
    
    # Filtering to get resistance data for the selected country of interest -> ggvis
    data_drug_resistance <- reactive({
      if (is.null(selected_country)){
        d <- gen_phen_country %>%
        group_by(country) %>%
        group_by(strain) %>%
        distinct() %>%
        gather(drug, resistance, rinh:rlin) %>%
        filter(resistance != "") %>% 
        select(drug, resistance) %>%
        mutate(ResistanceType = ifelse(resistance == "r", "1 - Resistant", 
                                       ifelse(resistance == "s", "2 - Sensitive", 
                                              "3 - Intermediate"))) %>%
        mutate(drug = substring(drug, 2)) %>% # remove 'r' prefix from drug name
        data.frame()
      } else {
        d <- gen_phen_country %>%
        group_by(country) %>%
        filter(country == as.character(selected_country$country[[1]][1])) %>%
        group_by(strain) %>%
        distinct() %>%
        gather(drug, resistance, rinh:rlin) %>%
        filter(resistance != "") %>% 
        select(drug, resistance) %>%
        mutate(ResistanceType = ifelse(resistance == "r", "1 - Resistant", ifelse(resistance == "s", "2 - Sensitive", "3 - Intermediate"))) %>%
        mutate(drug = substring(drug, 2)) %>% # remove 'r' prefix from drug name
        data.frame()
      }
    })
    
    # Labels and colors for Drug Resistance Charts
    fill_domain <- c("1 - Resistant", 
                     "2 - Sensitive", 
                     "3 - Intermediate")
    fill_range <- c(color_resistant_strains(), 
                    color_sensitive_strains(), 
                    color_intermediate_strains())
    
    # Hover Content for Drug Resistance Charts
    hover_drug_resistance <- function(x) { 
      if(is.null(x)){
        return(NULL)
      }
      paste0("<strong>", c("Resistance Type", "Count"), ": </strong>", 
             format(x)[c(1, 4)], collapse = "<br/>") 
    }
    
    data_drug_resistance %>%
      ggvis(~drug) %>%
      layer_bars(fill = ~ResistanceType, opacity := .85) %>%
      scale_ordinal("fill", domain = fill_domain, range = fill_range) %>%
      add_axis("x", title = "",
               properties = axis_props(labels = list(angle = 45,
                                                     align = "left")
               )) %>% 
      add_tooltip(hover_drug_resistance, "hover") %>% 
      bind_shiny("resistance_plot", "resistance_plot_ui")
    
    data_drug_resistance_gene <- reactive({
      if (is.null(selected_country)){
        summary_gen_phen_country
      } else {
        subset(data.frame(summary_gen_phen_country),
               country == as.character(selected_country$country[[1]][1]))[, !colnames(summary_gen_phen_country) %in% c("code", "country")]
      }
    })
    
    output$country1 <- renderText({
      as.character(selected_country$country[[1]][1])
    })
    
    output$country2 <- renderText({
      as.character(selected_country$country[[1]][1])
    })
    
    output$country3 <- renderText({
      as.character(selected_country$country[[1]][1])
    })

    output$download_country_data <- downloadHandler(
      filename = function() { paste("tb_resistance_",
                                    as.character(selected_country$country[[1]][1]),
                                    '.txt', sep = '') },
      content = function(file) {
        write.csv(data_drug_resistance(), file)
      }
    )
    
    output$download_world_data <- downloadHandler(
      filename = function() { paste("tb_resistance_world",
                                    '.txt', sep = '') },
      content = function(file) {
        write.csv(data_drug_resistance(), file)
      }
    )
    
    selected_country_lineage <- reactive({ 
      country_name_holder <- NULL
      if (!is.null(selected_country)) {
        country_name_holder <- selected_country
        selected_country_lineage <- subset(gen_phen_country_by_strain, 
                                           country == as.character(country_name_holder$country[[1]][1]))
      } else {
        country_name_holder <- NULL
        selected_country_lineage <- gen_phen_country_by_strain
      }
      data.frame(selected_country_lineage)
      
    })
    
    # Hover Content for Lineage Charts
    hover_lineage <- function(x) { 
      if(is.null(x)){
        return(NULL)
      }
      paste0("<strong> Count: </strong>", format(x)[3], collapse = "<br/>") 
    }
    
      lineage_data <- subset(selected_country_lineage(), 
                           !is.na(family) & spfamily_parentstrain != "" & spfamily_parentstrain != "NA")
      
      if(nrow(lineage_data) == 0){
        output$lineage_indicator <- renderText({
          "  "
        })
        output$error_message <- renderText({
          "No Lineage Data Available"
        })
      } else {
        output$lineage_indicator <- renderText({
          " "
        })
        lineage_data %>%
          ggvis(~spfamily_parentstrain, fill := "lightgrey") %>%
          layer_bars() %>%
          add_axis("x", title = "",
                   properties = axis_props(labels = list(angle = 45,
                                                         align = "left")
                   )) %>%
          add_tooltip(hover_lineage, "hover") %>%
          bind_shiny("lineage_plot", "lineage_plot_ui")
        
        output$error_message <- renderText({
          ""
        })
      }
    
    # Populating user choices for Mutations charts
    observe({
      updateSelectizeInput(session, 'selected_locus',
                           'Select a locus', 
                           choices = sort(unique(subset(data.frame(summary_gen_phen_country),
                                                        country == as.character(selected_country$country[[1]][1]))$locus)),
                           options = list(
                             placeholder = 'Select a locus',
                             onInitialize = I('function() { this.setValue(""); }')))
    })

# @asiegmann: Need a Ninja to figure out how to make this work, and then implement the actual filtering in mutation_function
    observe({
        updateSelectizeInput(session, 'selected_variant',
                             'Select a variant', 
                             choices = sort(unique(subset(data.frame(summary_gen_phen_country),
                                                          country == as.character(selected_country$country[[1]][1]))$variant)),
                             options = list(
                                 placeholder = 'Select a variant',
                                 onInitialize = I('function() { this.setValue(""); }')))
    })
    
    observe({
      updateSelectizeInput(session, 'selected_drug',
                           'Select a drug',
                           choices = sort(c("inh - Isoniazid" = "rinh",    
                                       "eth - Ethionamide" = "reth",
                                       "rif - Rifampicin"= "rrif",
                                       "rfb - Rifabutin" = "rrfb",   
                                       "emb - Ethambutol" = "remb",
                                       "pza - Pyrazinamide"= "rpza",
                                       "str - Streptomycin" = "rstr",
                                       "amk - Amikacin" = "ramk",
                                       "cap - Capreomycin" = "rcap",
                                       "kan - Kanamycin" = "rkan",
                                       "oflx - Ofloxacin" = "roflx",
                                       "cip - Ciprofloxacin" = "rcip",
                                       "levo - Levofloxacin" = "rlevo",
                                       "pas - Para-Aminosalicylic Acid" = "rpas",
                                       "cys - Cycloserine" = "rcys",
                                       "tha - Thioacetazone" = "rtha",
                                       "pro - Prothionamide" = "rpro",
                                       "clof - Clofazimine" = "rclof",
                                       "moxi - Moxifloxacin" = "rmoxi",
                                       "clar - Clarithromycin" = "rclar",
                                       "gati - Gatifloxacin" = "rgati",
                                       "amoxclav - Amoxicillin Clavulanate" = "ramoxclav",
                                       "lin - Linezolid" = "rlin")),
                           options = list(
                             placeholder = 'Select a drug',
                             onInitialize = I('function() { this.setValue(""); }')))
    })
    
    # Fetching the country-specific mutation data
    selected_country_mutations <- reactive({ 
      country_name_holder <- NULL
      if (!is.null(selected_country)) {
        country_name_holder <- selected_country
        selected_country_mutations <- subset(tidy_gen_phen_country, 
                                             country == as.character(country_name_holder$country[[1]][1]))
      }
      else {
        country_name_holder <- NULL
        selected_country_mutations <- gen_phen_country_by_strain
      }
      data.frame(selected_country_mutations)
    })
    
    fill_domain <- c("r", 
                     "s", 
                     "i") 
    
    # Hover Content for Mutations Charts
    hover_mutations <- function(x) { 
      if(is.null(x)){
        return(NULL)
      }
      paste0("<strong>", c("Resistance Type", "Count"), ": </strong>", 
             format(x)[c(1, 4)], collapse = "<br/>") 
    }
    
    
    # Set default message for pre-select panel
    output$mutations_indicator <- renderText({
      "  "
    })
    output$error_message2 <- renderText({
      "Select a locus, variant, and drug above - then press 'GO'"
    })
    
    
    # Function to trigger when user presses 'GO' action button
    mutation_func <- function(mutations_data, selected_locus, selected_drug) {
      
        # Filtering selected_country_mutations by user inputted locus/variant/drug + removing blanks/NAs
        mutations_data <- selected_country_mutations() %>% 
          filter(locus == selected_locus)
        
        if (nrow(mutations_data) == 0) {
          output$mutations_indicator <- renderText({
            "  "
          })
          
          progress <- shiny::Progress$new()
          progress$set(message = "Fetching data", value = 0)
          # Close the progress when this reactive exits (even if there's an error)
          on.exit(progress$close())
          
          output$error_message2 <- renderText({
            "No locus-specific mutation data - please choose another combination and press 'GO'"
          })
        } else {
          mutations_data <- mutations_data %>%  filter(drug == selected_drug)
          
          if (nrow(mutations_data) == 0) {
            output$mutations_indicator <- renderText({
              "  "
            })
            
            progress <- shiny::Progress$new()
            progress$set(message = "Fetching data", value = 0)
            # Close the progress when this reactive exits (even if there's an error)
            on.exit(progress$close())
            
            output$error_message2 <- renderText({
              "No drug-specific mutation data - please choose another combination and press 'GO'"
            })
          } else {
            mutations_data <- mutations_data %>%
                              data.frame() %>% 
                              subset(!is.na(resistance) & resistance != "" & resistance != "NA")
            if (nrow(mutations_data) == 0) {
              output$mutations_indicator <- renderText({
                "  "
              })
              
              progress <- shiny::Progress$new()
              progress$set(message = "Fetching data", value = 0)
              # Close the progress when this reactive exits (even if there's an error)
              on.exit(progress$close())
              
              output$error_message2 <- renderText({
                  "No drug specific mutation data - please choose another combination and press 'GO'"
              })
            } else {
              output$mutations_indicator <- renderText({
                " "
              })
              mutations_data %>% 
                ggvis(~resistance) %>% 
                layer_bars(fill = ~resistance) %>%  
                scale_ordinal("fill", domain = fill_domain, range = fill_range) %>% 
                add_axis("x", title = "Resistance Type") %>%
                add_tooltip(hover_mutations, "hover") %>% 
                bind_shiny("mutations_plot", "mutations_plot_ui")

              progress <- shiny::Progress$new()
              progress$set(message = "Fetching data", value = 0)
              # Close the progress when this reactive exits (even if there's an error)
              on.exit(progress$close())
              
              output$error_message2 <- renderText({
                ""
              })
            }
          }
        }
    }
    
    observeEvent(input$mutation_button, {
      mutation_func(mutations_data, input$selected_locus, input$selected_drug)
      })
    
    })
  
})
