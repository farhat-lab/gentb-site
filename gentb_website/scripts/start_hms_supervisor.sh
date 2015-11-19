#!/bin/bash

# Set up the virtualenv
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
source venv_tb/bin/postactivate

cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website

# Start supervisor (not happening yet....)
supervisord -c /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings/hms_supervisord.conf

# Stop supervisor
# supervisorctl -c /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings/hms_supervisord.conf

# For crontab
#
#@reboot /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/scripts/start_hms_supervisor.sh
 
