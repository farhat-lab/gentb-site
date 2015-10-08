#!/bin/sh
# This script should be wrapped by another script that
# encloses all of these commands in "scl enable python27"
# and is run by the "plaid" user one time for setup.
# su plaid -l -s /bin/sh -c 'scl enable python27 "path/to/script.sh"'
# See also http://developerblog.redhat.com/2013/02/14/setting-up-django-and-python-2-7-on-red-hat-enterprise-6-the-easy-way/

#source /usr/bin/virtualenvwrapper.sh
source /opt/rh/python27/root/usr/bin/virtualenvwrapper.sh
#
# Setup virtualenv
echo "Setup virtualenv"
mkdir -p /webapps/virtualenvs
export WORKON_HOME=/webapps/virtualenvs
mkvirtualenv tb_website
workon tb_website
# Install requirements (pip)
echo "Install requirements (pip)"
cd /webapps/code/PhthisisRavens/phthisis_website
pip install -r requirements/production.txt
#
# Validate settings file
#
echo "Validate settings file"
cd /webapps/code/PhthisisRavens/phthisis_website/tb_website
python manage.py validate --settings=tb_website.settings.production
#
# Create sqlite database + initial tables
#
echo "Create sqlite database + initial tables"
python manage.py syncdb --noinput --settings=tb_website.settings.production
#
# Create directory for sqlite db
#
echo "Create directory for sqlite db"
mkdir -p /webapps/data/tb
mkdir -p /webapps/data/tb/sqlite

echo "Create www directories (media/static/wsgi-related)"
mkdir /var/www/tb/media # user uploads
mkdir /var/www/tb/static # images, js, css, etc.
mkdir /var/www/tb/tb # wsgi.py
cp /webapps/code/PhthisisRavens/phthisis_website/tb_website/tb_website/vagrant-centos-wsgi.py /var/www/tb/tb/wsgi.py
python manage.py collectstatic --noinput --settings=tb_website.settings.production
