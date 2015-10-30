#! /bin/bash
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate

# should set the export DJANGO_SETTINGS_MODULE
source venv_tb/bin/postactivate

cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/apps/dropbox_helper/

# Download files from newly loaded dropbox links
python dropbox_retrieval_runner.py
