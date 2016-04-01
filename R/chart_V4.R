library(jsonlite)
library(tidyr)
library(rCharts)  
#to install rCharts: 
#1: install.packages('devtools') 
#2: library(devtools)
#3: install_github("ramnathv/rCharts")
library(ggvis)


##----- Select strain and drug

#selected_drug <- "inh"
#selected_strain <- "01-R1466"
#file_json<-'matrix.json'

create_scatterplot<-function(file_json) {
  ##----- Drug list
  drug_list <- c('inh','rif','pza', 'emb','str','eth','kan', 'cap',
                 'amk', 'cip', 'levo', 'oflx', 'pas') 
  ##----- Genetic region by drug
  gr <- list() # genetic region
  gr$inh <- c("katG", "inhA-promoter", "embB", "inhA", "iniB",
              "kasA", "ahpC", "embAB-promoter", "fabG1", "ndh", "oxyR")
  gr$rif <- c("rpoB")
  gr$pza <- c("pncA", "pncA-promoter")
  gr$emb <- c("embB", "embA", "embC", "embAB-promoter", "iniB", "iniA", "iniC")
  gr$str <- c("rpsL", "rrs", "gid", "rrs-rrl-promoter")
  gr$eth <- c("inhA-promoter", "ethA", "inhA", "fabG1")
  gr$kan <- c("rrs", "tlyA", "rrs-rrl-promoter")
  gr$cap <- c("rrs", "tlyA", "rrs-rrl-promoter")
  gr$amk <- c("rrs", "tlyA", "rrs-rrl-promoter")
  gr$cip <- c("gyrA", "gyrB")
  gr$levo <- c("gyrA", "gyrB")
  gr$oflx <- c("gyrA", "gyrB")
  gr$pas <- c("thyA")
  ##----- Read predict output from JSON file
  l <- fromJSON(txt = readLines(file_json))
  ##----- Important mutations
  important_mutations <- l[[2]]
  ##----- Construct important mutation matrix
  matrix_important_mutations <- NULL
  for (i in 1:length(important_mutations)) {
    mat <- data.frame(important_mutations[[i]])
    colnames(mat) <- drug_list
    mat <- gather(mat, drug, mutation, inh:pas)
    mat <- mat[!is.na(mat$mutation), ]
    mat$strain <- names(important_mutations[i])
    for (i in 1:nrow(mat))
      mat$genetic_region[i] <- tail(strsplit(as.character(mat$mutation)[i], "_")[[1]], 1)
    matrix_important_mutations <- rbind(matrix_important_mutations, mat)
  } 
  #matrix_important_mutations
  ##----- OTHER
  ##----- Other mutations
  other_mutations <- l[[3]]
  ##----- Construct other mutation matrix
  matrix_other_mutations <- NULL
  for (i in 1:length(other_mutations)) {
    mat <- data.frame(other_mutations[[i]])
    colnames(mat) <- drug_list
    mat <- gather(mat, drug, mutation, inh:pas)
    # mat <- mat[!is.na(mat$mutation), ]
    mat$strain <- names(other_mutations[i])
    for (i in 1:nrow(mat))
      mat$genetic_region[i] <- tail(strsplit(as.character(mat$mutation)[i], "_")[[1]], 1)
    matrix_other_mutations <- rbind(matrix_other_mutations, mat)
  } 
  #matrix_other_mutations
  ##----- Loop on strain and drug
  for (selected_drug in drug_list) {
    for (selected_strain in names(important_mutations)) {
      ##----- Create a count matrix
      count_important <- data.frame(matrix(gr[[selected_drug]],length(gr[[selected_drug]])))
      colnames(count_important) <- "genetic_region"
      count_important$count <- 0
      count_important$tip <- 0
      ##----- Subset by strain and drug
      selected_important <- subset(matrix_important_mutations,strain == selected_strain & drug == selected_drug)
      for (i in 1:nrow(count_important)) {
        ss <- subset(selected_important, genetic_region == count_important$genetic_region[i])
        count_important$count[[i]] <- nrow(ss)
        count_important$tip[[i]] <- paste(ss$mutation, collapse = ", ")
      }
      count_important$drug <-  selected_drug
      count_important$strain <- selected_strain
      ##----- Subset by strain and drug
      selected_other <- subset(matrix_other_mutations,strain == selected_strain & drug == selected_drug)
      ##----- Create a count matrix
      count_other <- data.frame(matrix(gr[[selected_drug]],length(gr[[selected_drug]])))
      colnames(count_other) <- "genetic_region"
      count_other$count <- 0
      #count.fields_other$tip <- 0
      for (i in 1:nrow(count_other)) {
        ss <- subset(selected_other, genetic_region == count_other$genetic_region[i])
        count_other$count[[i]] <- nrow(ss)
        count_other$tip[[i]] <- paste(ss$mutation, collapse = ", ")
      }
      count_other$drug <-  selected_drug
      count_other$strain <- selected_strain
      ##----- Combine IMPORTANT and OTHER
      count_important$type <- "Important"
      count_other$type <- "Other"
      ALL <- rbind(count_important, count_other)
      ##----- Graph
      graph <- rPlot(count ~ genetic_region, data = ALL, type = "point", color = "type",
                   tooltip = "#!function(item){ return item.tip }!#")
      graph$guides(x = list(title = "Genetic region"))
      graph$guides(y = list(title = "Number of mutations"))
      graph$guides(y = list(min = 0, max = 5))
      graph$guides(x = list(numticks = length(gr[[selected_drug]])))
      graph$guides(x = list(levels = gr[[selected_drug]]))
      # graph
      ##----- Save graph to HTML file
      graph$save(paste("scatterplot/scatterplot", selected_strain, selected_drug, "html", sep = "."), standalone = TRUE)
    }
  }
}
supressWarnings(create_scatterplot('matrix.json'))
