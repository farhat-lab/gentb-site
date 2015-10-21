# Overview of set-up on the [Orchestra Web Hosting Service](https://wiki.med.harvard.edu/Orchestra/WebHosting)

The https://gentb.hms.harvard.edu/ ```docroot``` uses an .htaccess file to reroute requests to the genTB Django app running via gunicorn.

The directory structure is as follows:

1. ```/www/gentb.hms.harvard.edu/docroot/```
    - contains .htaccess file
1. ```/www/gentb.hms.harvard.edu/docroot/tb/static```
  - Static files for Django app
1. ```/www/gentb.hms.harvard.edu/code/```
  - Repository with Django app
  - Contains Virtualenv running python 2.7.x
  - Gunicorn is run from ```flexatone.orchestra:9001```

## Set-up steps

- Log ```ssh -l username orchestra.med.harvard.edu```
- Create these directories:
```
cd /www/gentb.hms.harvard.edu/
mkdir code
# On the docroot
#
mkdir docroot/tb
mkdir docroot/tb/static
mkdir docroot/tb/media
# Update perms
#
chmod 755 -R docroot/tb
```

### Pull down the github repository

- Make sure you have an ssh key set up on the HMS server: https://help.github.com/articles/generating-ssh-keys/

```
cd /www/gentb.hms.harvard.edu/code
git clone git@github.com:IQSS/gentb-site.git
```

### Set up virtualenv

fyi: HMS server has virtualenv but not virtualenvwrapper

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
#
# Load python 2.7
module load dev/python/2.7.6
#
# Create virtualenv
virtualenv venv_tb
source venv_tb/bin/activate
#
# Load requirements, includes django and gunicorn
pip install -r requirements.txt
```

mysql connector installed manually
```
pip install --allow-external mysql-connector-python mysql-connector-python
```


### Reuse the virtualenv -- after setup

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
```

### Add production settings

Add these settings files:
  - ```secret_settings_prod_hms.json```
To this directory:
  - ```/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings```

### Check the settings

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website
python manage.py check --settings=tb_website.settings.production_hms
#
# Updates to prod settings
vim tb_website/settings/secret_settings_prod_hms.json
```

If this check fails, it is likely a setting in the ```secret_settings_prod_hms.json``` referred to above

### Create db tables and superuser

- Create tables
```
python manage.py migrate --settings=tb_website.settings.production_hms
```

- Create superuser
```
python manage.py createsuperuser --settings=tb_website.settings.production_hms
```

- Collect static files
This moves css, js, images, etc to the docroot.  Example image:
 - https://gentb.hms.harvard.edu/tb/static/images/TwoRavens.png

```
python manage.py collectstatic --settings=tb_website.settings.production_hms
```
Should see something like this:
  - ```109 static files copied to '/www/gentb.hms.harvard.edu/docroot/tb/static'.```
  - the 109 may be a different number


### Load Explore fixtures

These are the links to the Shiny server.

```
python manage.py loaddata apps/explore/fixtures/initial_data.json --settings=tb_website.settings.production_hms
```

### Set up the .htaccess file

Create an .htaccess file
```
vim /www/gentb.hms.harvard.edu/docroot/.htaccess
```

Add the following content:
```
# ----------------------------------
# Redirect requests to the Django app running on flexatone
# Gunicorn is being used to serve Django on http://flexatone.orchestra:9001
# -----------------------------------
RewriteEngine On
RewriteBase /
# -------------------------------
#
# Serve static files (css, js, images) directly
# These files live under /docroot/tb/...
# --------------------------------
RewriteCond %{REQUEST_URI} !^/tb/
# -------------------------------
#
# Send other requests to the Django app on flexatone
# -------------------------------
#ReWriteRule ^(.*)$ http://flexatone.orchestra:9001/$1 [P]
ReWriteRule ^(.*)$ http://gentb-app-prod01.orchestra:9001/$1 [P]# -------------------------------
#
# Temp redirect, if needed
# -------------------------------
#ReWriteRule ^(.*)$ /tb/static/images/predict.png
```

### Email Settings

Email message are sent via Django's [send_mail function](https://docs.djangoproject.com/en/1.8/topics/email/#send-mail).

Currently the project uses a Gmail account with settings supplied in the file: ```secret_settings_prod_hms.json```

Note: In order to work currently, the Gmail account must be set up with 2-step verification and an application password. For more information, read [How to generate an App password](https://support.google.com/accounts/answer/185833?hl=en)

To run a command line email test, log into the server and:
```


```


HMS QUESTIONS

flexatone - mysqldb driver
ImproperlyConfigured: Error loading MySQLdb module: libssl.so.10: cannot open shared object file: No such file or directory



### Move st
