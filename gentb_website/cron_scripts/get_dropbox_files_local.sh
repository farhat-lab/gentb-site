#! /bin/bash    

# invoke virtualenv
export WORKON_HOME=/Users/rmp553/.virtualenvs
source /usr/local/bin/virtualenvwrapper.sh
workon gentb

# run the script
cd /Users/rmp553/Documents/iqss-git/gentb-site/gentb_website//tb_website/apps/dropbox_helper/
python dropbox_retrieval_runner.py

