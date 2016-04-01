# devtools::install_github('ramnathv/rCharts')
library(rCharts)
library(jsonlite)

n1 <- rPlot(mpg ~ wt, data = mtcars, color = "gear", type = "point")

data(iris)
names(iris) = gsub('\\.', '', names(iris))
rPlot(SepalLength ~ SepalWidth | Species, data = iris, type = 'point', color = 'Species')

n1

n1$addControls("x", value = "wt", values = names(mtcars))
n1$addControls("y", value = "wt", values = names(mtcars))
n1$addControls("color", value = "gear", values = names(mtcars))
n1

n1$save("chart.html")

druglist <- c('inh','rif','pza', 'emb','str','eth','kan', 'cap', 'amk', 'cip', 'levo', 'oflx', 'pas') 

genetic_region <- c("katG", "inhA_promoter", "embB", "inhA"", iniB",
                    "kasA", "ahpC", "embAB_promoter", "fabG1", "ndh", "oxyR")

data_set <- data.frame(matrix(0, length(genetic_region), 4))
rownames(data_set) <- genetic_region
colnames(data_set) <- c("important",  "important_tip", "other", "other_tip")

data_set[1, 1] <- 2
data_set[1, 2] <- paste("SNP_CN_2155168_C944G_S315T_katG", "SNP_CN_4247431_G918A_M306I_embB")
data_set$genetic_region <- genetic_region

n1 <- rPlot(important ~ genetic_region, data = data_set, type = "point")
n1
n1$guides(x = list(levels =genetic_region))
n1

n1$save("chart.html")


p1 <- nPlot(mpg ~ wt, group = 'cyl', data = mtcars, type = 'multiBarChart')
p1$xAxis(axisLabel = 'Weight (in lb)')
p1

p2a <- nPlot(Freq ~ Hair, group = 'Eye', 
             data = hair_eye,
             type = 'multiBarChart'
)
p2a$chart(color = c('brown', 'blue', '#594c26', 'green'))
p2a$addFilters("Sex")
p2a$set(dom = 'chart2', width = 600)
p2a

levels(data_set$genetic_region) <- data_set$genetic_region

require(rCharts)
data(tips, package = 'reshape2')
r1 <- rPlot(tip ~ day, data = tips, type = 'point')
r1
n1

# library(tidyr)

s <- "SNP_CN_2155168_C944G_S315T_katG"
tail(strsplit(s, "_")[[1]], 1)


s <- "SNP_CN_4408156_A47C_L16R_gid"

file_json = "matrix.json"
l <- fromJSON(txt = readLines(file_json))
## important
ll <- l[[2]]
ll
names(ll)

lll <- ll[[1]]
is.matrix(lll)
colnames(lll) <- druglist
mutations <- na.omit(data.frame(lll[, selected_drug]))
colnames(mutations) <- "mutations"

mutations$genetic_region <- NA

for (i in 1:nrow(mutations))
  mutations$genetic_region[i] <- tail(strsplit(as.character(mutations$mutations)[i], "_")[[1]], 1)

mutations

## other
l3 <- l[[3]]
l3
names(l3)

# for (i in 1:nrow(mutations))
#   mutations$resistant[i] <- mutations$genetic_region[i] %in% genetic_region

selected_drug <- "inh"
selected_strain <- "01-R1466"
lll <- ll[[selected_strain]]
colnames(lll) <- druglist
lll[, selected_drug]

data_set <- NULL
for (i in length(ll)) 
  data_set <- rbind(data_set, ll[[1]])

names(d) <- c("strain", "drug","probability of resistance")

