# notes on shinyapps, mac install

### Log in for shinyapps.io
```
https://www.shinyapps.io/admin/#/login
```

- **Note**: The application instance size should be: ```X-Large (2GB)```
  - With the default 1GB, the maps were *not* showing up


### Following these instructions:
```
http://shiny.rstudio.com/articles/shinyapps.html
```

- They include (in R console):
```R
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("rstudio/shinyapps")

#install.packages('devtools')
devtools::install_github('rstudio/rsconnect')
library(rsconnect)
```


### Package installs (for local)


- Deploy commands

```R
library(rsconnect)
deployApp()
```

### Install list for local runs:

```R
install.packages('shiny')
install.packages('leaflet')
install.packages('RColorBrewer')
install.packages('maps')
install.packages('data.table')
install.packages('dplyr')
install.packages('ggvis')
install.packages('tidyr')
install.packages('ggplot2')
install.packages('DT')
install.packages('stringr')
install.packages('countrycode')
```
