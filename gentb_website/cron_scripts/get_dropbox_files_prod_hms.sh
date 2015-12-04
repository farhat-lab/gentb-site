#! /bin/bash

# Activate the virtual environment
#
source /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/activate

# Should set the export DJANGO_SETTINGS_MODULE
#
source /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/postactivate

# Download files from newly loaded dropbox links
#
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/python /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/apps/dropbox_helper/dropbox_retrieval_runner.py
