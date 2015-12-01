library(shiny)
library(ggvis)
library(tidyr)

shinyServer(function(input, output) {
  
  output$distPlot <- renderPlot({
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
  })
})

# all_values <- function(x) {
#   if (is.null(x))
#     return(NULL)
#   paste0(format(x)[c(1, 4)])
# }
# 
# data_set
# data_set["embB", "other"] <- 3
# 
# data_set
# 
# d <- data_set %>%
#   gather(type, number, important, other)
# 
# d$tip <- paste(as.character(d$important_tip))
# d$id <- 1:nrow(d)
# 
# ##----- Get tooltip
# # http://stackoverflow.com/questions/24519980/add-data-to-ggvis-tooltip-thats-contained-in-the-input-dataset-but-not-directly/24528087#24528087
# 
# all_values <- function(x) {
#   if(is.null(x)) return(NULL)
#   row <- d[d$id == x$id, ]
#   # print(row)
#   row$id <- NULL
#   row$important_tip <- NULL
#   row$value <- NULL
#   row$other_tip <- NULL
#   paste0(names(row), ": ", format(row), collapse = "<br />")
# }
# 
# d %>%
#   ggvis(~genetic_region, ~number, fill = ~type, key := ~id) %>%
#   layer_points() %>%
#   add_tooltip(all_values, "hover")
# 
