#! /bin/bash

# invoke virtualenv
export WORKON_HOME=/Users/rmp553/.virtualenvs
source /usr/local/bin/virtualenvwrapper.sh
workon gentb

# Set the django module
#
export DJANGO_SETTINGS_MODULE=tb_website.settings.local

# Run the pipeline for the next available PredictDataset
#
cd /Users/rmp553/Documents/iqss-git/gentb-site/gentb_website/tb_website/apps/predict/
python pipeline_hardcoded_script_runner.py
