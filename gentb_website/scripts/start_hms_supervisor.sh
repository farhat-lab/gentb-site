#!/bin/bash

# Set up the virtualenv
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
source venv_tb/bin/postactivate

# Start supervisor (not happening yet....)
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/gunicorn -c tb_website/settings/hms_gunicorn_config.py tb_website.wsgi:application
