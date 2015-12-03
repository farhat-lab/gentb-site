Notes on shinyapps.io
=====================

## Links

### Viewing
  - Current url: https://hmdc.shinyapps.io/genTB
  - genTB iframed view: https://gentb.hms.harvard.edu/maps/tb-map/

### Login/Config in for shinyapps.io

  - https://www.shinyapps.io/admin/#/login
  - **Note**: Make sure the application instance size is: ```X-Large (2GB)```
    - With the default 1GB, the maps were *not* showing up

## Update/Deploy

Note, this was done on a Mac using RStudio.
  - RStudio, version 0.99.486
  - R, version 3.2.2

### Official RStudio instructions:
```
http://shiny.rstudio.com/articles/shinyapps.html
```

### Initial Setup

1. Start RStudio
1. In the R console, run preliminaries from instructions above.  
  - This includes installing RStudio packages:
  ```R
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("rstudio/shinyapps")

#install.packages('devtools')
devtools::install_github('rstudio/rsconnect')
library(rsconnect)
```
3. Set working directory using this github repository
  - Menu: Session -> Set working directory -> Choose directory
    - choose ```(your local dir)/gentb-site/shinyapps.io/genTB```


### Install packages for local runs

Run these commands via the R console:

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

### Run/Test it locally

  - Open the "ui.R" file
    - location: ```(your local dir)/gentb-site/shinyapps.io/genTB/ui.R```
  - Click "Run App" (top right of file editor)
  - Does it work?

### Deploy to server!

  - Run these two commands in the R studio console

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

- View the output: https://hmdc.shinyapps.io/genTB
