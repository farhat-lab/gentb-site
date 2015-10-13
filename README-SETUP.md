# Overview of set-up on the [Orchestra Web Hosting Service](https://wiki.med.harvard.edu/Orchestra/WebHosting)

The https://gentb.hms.harvard.edu/ ```docroot``` uses an .htaccess file to reroute requests to the genTB Django app running via gunicorn.

The directory structure is as follows:

1. ```/www/gentb.hms.harvard.edu/docroot/```
    - contains .htaccess file
1. ```/www/gentb.hms.harvard.edu/docroot/static-files```
  - Static files for Django app
1. ```/www/gentb.hms.harvard.edu/gentb-app/```
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

```
cd /www/gentb.hms.harvard.edu/code/gentb-site/gentb_website
# Load python 2.7
module load dev/python/2.7.6
# Create virtualenv
virtualenv venv_tb
source venv_tb/bin/activate
# Load requirements
pip install -r requirements.txt
```
