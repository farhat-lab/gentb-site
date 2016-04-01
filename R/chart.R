library(jsonlite)
library(tidyr)
library(rCharts)

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

file_json <- "matrix.json"
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

matrix_important_mutations

##----- Select strain and drug

IMPORTANT <- NULL

selected_drug <- "inh"
selected_strain <- "01-R1466"

for (selected_strain in names(important_mutations)) {
  
  for (selected_drug in drug_list) {
    
    ##----- Subset by strain and drug
    
    selected_important <- subset(matrix_important_mutations,
                                 strain == selected_strain & drug == selected_drug)
    
    ##----- Create a count matrix
    
    count_important <- data.frame(matrix(gr[[selected_drug]],
                                         length(gr[[selected_drug]])))
    colnames(count_important) <- "genetic_region"
    
    count_important$count <- 0
    count_important$tip <- 0
    
    for (i in 1:nrow(count_important)) {
      ss <- subset(selected_important, genetic_region == count_important$genetic_region[i])
      count_important$count[[i]] <- nrow(ss)
      count_important$tip[[i]] <- paste(ss$mutation, collapse = ", ")
    }
    count_important$drug <-  selected_drug
    count_important$strain <- selected_strain
    
    IMPORTANT <- rbind(IMPORTANT, count_important)
  }
}

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

matrix_other_mutations

##----- Select strain and drug

OTHER <- NULL

for (selected_strain in names(important_mutations)) {
  
  for (selected_drug in drug_list) {
    
    ##----- Subset by strain and drug
    
    selected_other <- subset(matrix_other_mutations,
                             strain == selected_strain & drug == selected_drug)
    
    ##----- Create a count matrix
    
    count_other <- data.frame(matrix(gr[[selected_drug]],
                                     length(gr[[selected_drug]])))
    colnames(count_other) <- "genetic_region"
    
    count_other$count <- 0
    count_other$tip <- 0
    
    for (i in 1:nrow(count_other)) {
      ss <- subset(selected_other, genetic_region == count_other$genetic_region[i])
      count_other$count[[i]] <- nrow(ss)
      count_other$tip[[i]] <- paste(ss$mutation, collapse = ", ")
    }
    count_other$drug <-  selected_drug
    count_other$strain <- selected_strain
    
    OTHER <- rbind(OTHER, count_other)
  }
}

##----- Combine IMPORTANT and OTHER

IMPORTANT$type <- "Important"
OTHER$type <- "Other"

ALL <- rbind(IMPORTANT, OTHER)

ALL

##----- Graph

library(ggvis)

d <- subset(ALL, drug == "inh" 
            & genetic_region %in% gr[["inh"]] 
            & strain == selected_strain)
d <- droplevels(d)
levels(d$genetic_region) <- gr$inh

d %>%
  ggvis(~genetic_region, ~count) %>%
  layer_points()

# n1 <- rPlot(count ~ genetic_region, data = IMPORTANT, type = "point")
# n1
# selected_strain <- "01-R1466"
# selected_drug <- "inh"

n1 <- rPlot(count ~ genetic_region, data = droplevels(IMPORTANT), type = "point")
n1

n1$addControls("drug", value = "wt", values = drug_list)

n1
n1$addControls("y", value = genetic_region, values = names(mtcars))
n1$addControls("color", value = "gear", values = names(mtcars))

iris %>% 
  ggvis(x = input_select(c('Petal.Width', 'Sepal.Length'), map = as.name)) %>% 
  layer_points(y = ~Petal.Length, fill = ~Species)
