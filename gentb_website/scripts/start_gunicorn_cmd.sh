#!/bin/bash
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
source venv_tb/bin/postactivate
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/gunicorn -c tb_website/settings/hms_gunicorn_config.py tb_website.wsgi:application
