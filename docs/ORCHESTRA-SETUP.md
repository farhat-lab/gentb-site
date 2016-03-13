([go to code update section](#code-updates))

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
  - Gunicorn is run from ```rc-app-shared01.orchestra:9001```

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
# created the "code" directory if it doesn't exist
cd /www/gentb.hms.harvard.edu/code
git clone git@github.com:IQSS/gentb-site.git
```

### Logging

#### Set up logging directories

```
mkdir /www/gentb.hms.harvard.edu/logging
mkdir /www/gentb.hms.harvard.edu/logging/gentb-django
mkdir /www/gentb.hms.harvard.edu/logging/supervisor
mkdir /www/gentb.hms.harvard.edu/logging/gunicorn
```

#### Where are the log settings?

1. *Django*:
    - Log settings: ```settings/production_hms.py```.  
        - At the time of this writing, logs are set to rotate, with the last 5 being kept.  Each log has a max size of 5mb.
    - Tailing the log:
        - ```tail /www/gentb.hms.harvard.edu/logging/gentb-django/gentb.log```
2. *Gunicorn*:
    - Log settings: ```settings/hms_gunicorn_config.py```
    - Tailing the logs:
        - Access: ```tail /www/gentb.hms.harvard.edu/logging/gunicorn/gunicorn_access.log```
        - Error: ```tail /www/gentb.hms.harvard.edu/logging/gunicorn/gunicorn_error.log```
3. *Supervisor*:
    - Log settings: ```settings/hms_supervisord.conf```
        - At the time of this writing, logs are set to rotate, with the last 5 being kept.  Each log has a max size of 10mb.
    - Tailing the log:
        - ```tail /www/gentb.hms.harvard.edu/logging/supervisor/supervisord.log```


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

### Make postactivate script

Note virtualenvwrapper currently isn't available on the orchestra system.  This is a DIY postactivate script

1.  Create a ```postactivate``` file
```
vim venv_tb/bin/postactivate
```
1. Add these lines to the file:
```
#!/bin/bash
# This hook is run after this virtualenv is activated.
export DJANGO_SETTINGS_MODULE=tb_website.settings.production_hms
. /opt/lsf/conf/profile.lsf
```
1.  Run the ```postactivate``` script
```
source venv_tb/bin/postactivate
```

### Reuse the virtualenv -- after setup

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
source venv_tb/bin/postactivate
# If profile not set to run bsub
 . /opt/lsf/conf/profile.lsf
```

### Add production settings

Add this settings file:
  - ```secret_settings_prod_hms.json```
To this directory:
  - ```/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings```
  - Full path: ```/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings/secret_settings_prod_hms.json```

### Check the settings

- Note: If you create/run the ```postactivate``` script above, then the ```--settings=...``` option is not needed

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

  - Note: If you create/run the ```postactivate``` script above, then the ```--settings=...``` option is not needed


- Create superuser
```
python manage.py createsuperuser --settings=tb_website.settings.production_hms
```

  - Note: If you create/run the ```postactivate``` script above, then the ```--settings=...``` option is not needed


- Collect static files
This moves css, js, images, etc to the docroot.  Example image:
 - https://gentb.hms.harvard.edu/tb/static/images/TwoRavens.png

```
python manage.py collectstatic --settings=tb_website.settings.production_hms
```
Should see something like this:
  - ```109 static files copied to '/www/gentb.hms.harvard.edu/docroot/tb/static'.```
  - the 109 may be a different number



### Set the site name

```
python manage.py shell --settings=tb_website.settings.production_hms
from django.contrib.sites.models import Site
s = Site.objects.all()[0]
s.name = 'gentb.hms.harvard.edu'
s.domain = 'gentb.hms.harvard.edu'
s.save()
```

### Load initial status settings + pipeline script directory

```
# status settings
python manage.py loaddata apps/predict/fixtures/predict_statuses.json

# location of pipeline scripts (HMS)
python manage.py loaddata apps/predict/fixtures/hms_pipeline_scripts_dir.json
```



### Set up the .htaccess file

Create an .htaccess file
```
vim /www/gentb.hms.harvard.edu/docroot/.htaccess
```

Add the following content:
```
# ----------------------------------
# Redirect requests to the Django app running on rc-app-shared01
# Gunicorn is being used to serve Django on http://rc-app-shared01.orchestra:9001
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
# Send other requests to the Django app on rc-app-shared01
# -------------------------------
#ReWriteRule ^(.*)$ http://flexatone.orchestra:9001/$1 [P]
#ReWriteRule ^(.*)$ http://gentb-app-prod01.orchestra:9001/$1 [P]
ReWriteRule ^(.*)$ http://rc-app-shared01.orchestra:9001/$1 [P]
# -------------------------------
#
# Temp redirect, if needed
# -------------------------------
#ReWriteRule ^(.*)$ /tb/static/images/predict.png
```


### Email Settings

Email message are sent via Django's [send_mail function](https://docs.djangoproject.com/en/1.8/topics/email/#send-mail).

Currently the project uses a Gmail account with settings supplied in the file:
  - file: ```secret_settings_prod_hms.json```
  - prod location: ```/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings/secret_settings_prod_hms.json```

Note: In order to work currently, the Gmail account must be set up with 2-step verification and an application password. For more information, read [How to generate an App password](https://support.google.com/accounts/answer/185833?hl=en)

To run a command line email test, log into the server and go into the Django shell:

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
source venv_tb/bin/postactivate
cd tb_website
python manage.py shell
```

Next, from the django shell, try the ```send_mail``` command:

```
# !! put your email address in the line below
to_email = "--YOUR EMAIL ADDRESS HERE--"
#
from django.conf import settings
from django.core.mail import send_mail
send_mail(subject='genTB test email',
        message='test the email settings',
        from_email=settings.DEFAULT_FROM_EMAIL,
        recipient_list=[to_email])

```

If it works, you should receive an email shortly.  An error message will appear on failure.

### Run gunicorn/supervisord


#### Run supervisord
```
# Run a script that kicks off supervisord
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/scripts/start_hms_supervisor.sh
```

The ```start_hms_supervisor.sh``` script above contains these commands:
```
# Set up the virtualenv
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
source venv_tb/bin/postactivate

cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website

# Start supervisor
supervisord -c /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings/hms_supervisord.conf

# Start the supervisor control prompt
# e.g. issue 'status', 'shutdown', 'help', etc.
supervisorctl -c /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings/hms_supervisord.conf
```

##### Stop and Restart supervisord

-  **Stop it**
```
# Assumes in virtualenv above

cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website

# Activate the supervisor command line control
supervisorctl -c /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/tb_website/settings/hms_supervisord.conf
shutdown
```

-  **Start it**
```
# Run a script that kicks off supervisord
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/scripts/start_hms_supervisor.sh
```

#### Running gunicorn without supervisord
```
# start the virtualenv
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
source venv_tb/bin/postactivate
cd tb_website
gunicorn -c tb_website/settings/hms_gunicorn_config.py tb_website.wsgi:application

#gunicorn --log-file=- --workers=1 -b 0.0.0.0:9001 tb_website.wsgi:application &

```

- view processes (if needed)
```
ps -f -u <username>
```

### Load Explore fixtures (Skip, this has been moved to migrations)

These are the links to the Two Ravens app server.

```
python manage.py loaddata apps/explore/fixtures/initial_data.json --settings=tb_website.settings.production_hms
```

  - Note: If you create/run the ```postactivate``` script above, then the ```--settings=...``` option is not needed

### Crontabs

- Edit the crontab: ```crontab -e```
    - If on ```cron.orchestra```, use ```env EDITOR=vim crontab -e```
- Add the lines below, also found in the file ```gentb_website/cron_scripts/crontab```

```
# ------------------------------------------------
# Download files from new dropbox links (every 10 minutes)
# ------------------------------------------------
*/10 * * * * /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/cron_scripts/get_dropbox_files_prod_hms.sh
#
# ------------------------------------------------
# See if there are files ready for pipeline processing (every 15 minutes)
# ------------------------------------------------
*/15 * * * * /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/cron_scripts/run_pipeline_prod_hms.sh
#
# --------------------------------------
# On reboot, get supervisor running the genTB site again
# --------------------------------------
@reboot /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/scripts/start_hms_supervisor.sh
#
```

### Test script Run

Test ```ScriptToRun``` object to add via the admin:
```
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/python  /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/apps/script_helper/test_script.py
```

#### Example of remote test of callback:
    - called HMS serve from local machine
```
>>> url = 'https://gentb.hms.harvard.edu/predict/pipeline-run-results-notice/'
>>> p = {'result_data': '{"amk":{"r":[0.7683],"s":[0.2317]},"cap":{"r":[0.5202],"s":[0.4798]},"cip":{"r":[0.6477],"s":[0.3523]},"emb":{"r":[0.09],"s":[0.91]},"gai":{"r":[0.2632],"s":[0.7368]},"inh":{"r":[0.8933],"s":[0.1067]},"kan":{"r":[0.9487],"s":[0.0513]},"levo":{"r":[0.0042],"s":[0.9958]},"moxi":{"r":[0.3219],"s":[0.6781]},"oflx":{"r":[0.0091],"s":[0.9909]},"pas":{"r":[0.5854],"s":[0.4146]},"pza":{"r":[0.163],"s":[0.837]},"rif":{"r":[0.7605],"s":[0.2395]},"str":{"r":[0.8175],"s":[0.1825]}}', 'success': True, 'run_md5': u'2af2023f8c6103dbe00803929de54c28'}
>>> import requests
>>> requests.post(url, data=p)
<Response [200]>
>>>
```

#### check this:

```
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/venv_tb/bin/python  /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website/apps/script_helper/test_script.py '{"file_directory": "/home/rp188/gentb_test/tbdata_00000001", "run_md5": "2af2023f8c6103dbe00803929de54c28", "admin_url": "https://gentb.hms.harvard.edu/gentb-admin/predict/predictdataset/1/", "callback_url": "https://gentb.hms.harvard.edu/predict/pipeline-run-results-notice/", "dataset_id": 1, "user_email": "raman_prasad@harvard.edu"}'
```

# Code updates

## Login
```
ssh username@orchestra.med.harvard.edu
ssh username@rc-app-shared01.orchestra
```

## Reuse the virtualenv -- after setup

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
source venv_tb/bin/activate
source venv_tb/bin/postactivate
# If profile not set to run bsub
 . /opt/lsf/conf/profile.lsf
```

## Pull from github

```
git pull
```

## Update static files

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website
python manage.py collectstatic
```

## Run database schema updates

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website
python manage.py migrate
```

## Restart server

- See [Running supervisord](#run-gunicornsupervisord)

## Run dropbox retrieval script

```
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/cron_scripts/get_dropbox_files_prod_hms.sh
```

## Run pipeline script

```
/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/cron_scripts/run_pipeline_prod_hms.sh
```
