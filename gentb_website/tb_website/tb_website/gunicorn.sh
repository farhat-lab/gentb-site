#!/bin/bash
NAME="gentb_hms" # Name of the application
DJANGODIR=/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website # Django project directory
SOCKFILE=/home/osboxes/projects/active/django_project/gunicorn.sock # we will communicte using this unix socket

USER=rp188 # the user to run as
GROUP=gentb # the group to run as
NUM_WORKERS=3 # how many worker processes should Gunicorn spawn

MAX_REQUESTS=1 # reload the application server for each request
DJANGO_SETTINGS_MODULE=tb_website.settings.production_hms # which settings file should Django use
DJANGO_WSGI_MODULE=tb_website.wsgi:application # WSGI module name

echo “Starting $NAME as `whoami`"

# Activate the virtual environment
cd $DJANGODIR
source /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/activate

export DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE
export PYTHONPATH=$DJANGODIR:$PYTHONPATH

# Create the run directory if it doesn’t exist
RUNDIR=$(dirname $SOCKFILE)
test -d $RUNDIR || mkdir -p $RUNDIR

# Start your Django Unicorn

# Programs meant to be run under supervisor should not daemonize themselves (do not use –daemon)
exec /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/gunicorn ${DJANGO_WSGI_MODULE}:application \
–name $NAME \
–workers $NUM_WORKERS \
–max-requests $MAX_REQUESTS \
–user=$USER –group=$GROUP \
–bind=0.0.0.0:9000 \
–log-level=error \
–log-file=-
