#! /bin/bash
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate

# Should set the DJANGO_SETTINGS_MODULE
#
source venv_tb/bin/postactivate

# Run the pipeline for the next available PredictDataset
#
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/apps/predict/
python pipeline_hardcoded_script_runner.py --next
