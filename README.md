# genTB website

Basic django site to facilitate the contribution and analysis of fastQ and VCF files.

Docs:
 - Basic informational web pages (this README.md)
 - Set-up/Updates on Orchestra: [docs/README-SETUP-ORCHESTRA](docs/README-SETUP-ORCHESTRA.md)
 - Shiny Maps (shinyapp.io): [Set-up/Updates](shinyapps.io/NOTES-SHINY.md)
 - Predict workflow: [docs/README-PREDICT-WORKFLOW.md](docs/README-PREDICT-WORKFLOW.md)

---
## Authentication

Currently users the basic Django auth.

## Pages

### Homepage

A basic view/template combination

  - View: "view_homepage" [tb_website/apps/basic_pages/views.py](tb_website/apps/basic_pages/views.py)
  - Template: [tb_website/templates/homepage.html](tb_website/templates/homepage.html)

![homepage screenshot](screen-shots/genTB-home.png?raw=true "genTB Homepage")

### Predict

This pages requires the user to be logged in.  

### Data Page

A basic view/template combination.  The template includes embedded [Dataverse Widgets](http://datascience.iq.harvard.edu/blog/dataverse-40-theme-widgets) in the form of javascript snippets

  - View: "view_data_page" [tb_website/apps/basic_pages/views.py](tb_website/apps/basic_pages/views.py)
  - Template: [tb_website/templates/data_page.html](tb_website/templates/data_page.html)

![Share Page screenshot](screen-shots/genTB-data.png?raw=true "genTB Data page")


### Map

A basic view/template combination.  The maps are actually hosted by a Shiny application and embedded in an iframe.  The url to the iframe is written directly into the template.

  - View: "map" [tb_website/apps/maps/views.py](tb_website/apps/maps/views.py)
  - Template: [tb_website/templates/maps/basic_map.html](tb_website/templates/maps/basic_map.html)

- **Note**: The Shiny app url is currently https://hmdc.shinyapps.io/genTB and part of the IQSS Shinyapps.io service
    - For more information, see [shinyapps.io/NOTES-SHINY.md](shinyapps.io/NOTES-SHINY.md)

![Map Page screenshot](screen-shots/genTB-map.png?raw=true "genTB Map page")


### Explore (TwoRavens)

[TwoRavens](https://github.com/IQSS/TwoRavens) is a system of interlocking statistical tools for data exploration, analysis, and meta-analysis.  This page uses TwoRavens via an iframe to analyze the genTB master data file.

TwoRavens needs two values which are supplied via the genTB database:

  1. **codebook_file_url** - This url links directly to PDF file on Dataverse and identified by file id number
    - Example ```codebook_file_url``: https://dataverse.harvard.edu/api/access/datafile/2694344
  2. **two_ravens_url** - This url contains the TwoRavens link and is used within an iframe.  
    - The url contains both a **Dataverse file id number** and an **API key**
        - This API key within the url is bad practice and will probably be changed on the Dataverse side
    - Example ```two_ravens_url```: https://rserve.dataverse.harvard.edu/dataexplore/gui.html?dfId=2693726&key=c54f07b7-5098-461c-adf3-a976c0d62f6e

These values may be supplied via the Django admin:

    - Models.py file for the values above:
        - [tb_website/apps/explore/models.py](tb_website/apps/explore/models.py)

![Explore Admin Page screenshot](screen-shots/genTB-explore-admin.png?raw=true "genTB Explore Admin page")


The view/template files may be found here:

  - View: "explore:home" [tb_website/apps/basic_pages/views.py](tb_website/apps/basic_pages/views.py)
  - Template: [tb_website/templates/explore.html](tb_website/templates/explore.html)

![Explore Page screenshot](screen-shots/genTB-explore.png?raw=true "genTB Explore page")

### About Page

A basic view/template combination

  - View: "about" [tb_website/apps/basic_pages/views.py](tb_website/apps/basic_pages/views.py)
  - Template: [tb_website/templates/about.html](tb_website/templates/about.html)

![about page screenshot](screen-shots/genTB-about.png?raw=true "genTB About Page")
