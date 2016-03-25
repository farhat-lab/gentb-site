#!/bin/bash

# Go to this directory (where the script is held)
cd "$( dirname "${BASH_SOURCE[0]}" )"
cd ..

# Set up the virtualenv
source pythonenv/bin/activate

# Start gunicorn
gunicorn -c "$PWD/scripts/gunicorn_config.py" tb_website.wsgi:application
 
