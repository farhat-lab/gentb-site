# genTB website

Basic django site to facilitate the contribution and analysis of fastQ and VCF files.

---
## Authentication

Currently users the basic Django auth.

## Pages

### Homepage

A basic view/template combination

  - View: "view_homepage" [gentb_website/tb_website/apps/basic_pages/views.py](gentb_website/tb_website/apps/basic_pages/views.py)
  - Template: [gentb-site/gentb_website/tb_website/templates/homepage.html](gentb-site/gentb_website/tb_website/templates/homepage.html)

![TB homepage screenshot](screen-shots/genTB-home.png?raw=true "genTB Homepage")

### Share

A basic view/template combination.  The template includes embedded [Dataverse Widgets](http://datascience.iq.harvard.edu/blog/dataverse-40-theme-widgets) in the form of javascript snippets

  - View: "view_share_page" [gentb_website/tb_website/apps/basic_pages/views.py](gentb_website/tb_website/apps/basic_pages/views.py)
  - Template: [gentb-site/gentb_website/tb_website/templates/share.html](gentb-site/gentb_website/tb_website/templates/share.html)
