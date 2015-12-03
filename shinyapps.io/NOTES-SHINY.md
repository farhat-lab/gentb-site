# notes on shinyapps.io

- current url: https://hmdc.shinyapps.io/genTB

### Log in for shinyapps.io
```
https://www.shinyapps.io/admin/#/login
```

- **Note**: The application instance size should be: ```X-Large (2GB)```
  - With the default 1GB, the maps were *not* showing up


### Followed these instructions:
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

## RStudio on the Mac

### Set working directory
    - Start R Studio
    - Menu: Session -> Set working directory -> Choose directory
        - choose "(path to)/gentb-site/shinyapps.io/genTB"
### Install packages
    - Run the "Install list for local runs" above

### Run/Test it locally
    - Open the "ui.R" file under "...shinyapps.io/genTB"
    - Click "Run App" (top right of file editor)
    - Does it work?

### Deploy (same as above)

- Run these in the R studio console

```R
library(rsconnect)
deployApp()
```

- Example output:
    - Note: The warning message at the end is ok

```R
> deployApp()
> deployApp()
Preparing to deploy application...DONE
Uploading bundle for application: 65318...DONE
Deploying bundle: 323997 for application: 65318 ...
Waiting for task: 111256348
  building: Parsing manifest
  building: Building image: 318328
  building: Fetching packages
  building: Installing packages
  building: Installing files
  building: Pushing image: 318328
  deploying: Starting instances
  rollforward: Activating new instances
  terminating: Stopping old instances
Application successfully deployed to https://hmdc.shinyapps.io/genTB
Warning message:
invalid uid value replaced by that for user 'nobody'
```
