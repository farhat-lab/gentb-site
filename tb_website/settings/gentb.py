
import os
import sys
import json
import getpass

from tb_website.settings.base import *

SECRET_KEY = 'Nothing'

TB_SHARED_DATAFILE_DIRECTORY = ''

CONF_FILE  = os.path.expanduser(os.path.join('~', '.config', 'gentb-db.conf'))

def ask_for(name, slug, default=None, password=False):
    """Return a configuration item, or ask for new value"""
    try:
        return _ask(name, slug, default=default, password=password)
    except Exception as err:
        print "FAILED TO GET CONFIG ITEM: {} -> {}".format(name, err)

def _ask(name, slug, default=None, password=False):
    conf = {}
    if os.path.isfile(CONF_FILE):
        with open(CONF_FILE, 'r') as fhl:
            conf = json.loads(fhl.read())

    if '--new' in sys.argv:
        conf = {}

    if slug not in conf:
        if default is not  None:
            prompt = "{} ({}): ".format(name, default)
        else:
            prompt = "{}: ".format(name)

        if not password:
            ret = raw_input(prompt)
        else:
            ret = getpass.getpass(prompt)

        if not ret:
            if default is None:
                print("CANCELED on CONFIG ITEM {}".format(name))
                sys.exit(2)
            ret = default
        conf[slug] = ret

        CONFIG_DIR = os.path.dirname(CONF_FILE)
        if not os.path.isdir(CONFIG_DIR):
            os.makedirs(CONFIG_DIR)

        with open(CONF_FILE, 'w') as fhl:
            fhl.write(json.dumps(conf))

    return conf[slug]


DB_NAME = ask_for('MySQL DB Name or SQLite Filename', 'dbname')
DATABASES = {
  'default': {
     'ENGINE': 'django.db.backends.sqlite3',
     'NAME': DB_NAME,
  },
}
DB_IS_FILE = True

if '/' not in DB_NAME and not DB_NAME.endswith('.db'):
    DB_IS_FILE = False
    DATABASES['default']['ENGINE'] = 'django.db.backends.mysql'
    DATABASES['default']['HOST'] = ask_for('DB Host', 'dbhost', default='localhost')
    DATABASES['default']['USER'] = ask_for('DB User', 'dbuser')
    DATABASES['default']['PASSWORD'] = ask_for('DB Password', 'dbpass', password=True)

INSTALLED_APPS = list(set(INSTALLED_APPS) - set(WEBSITE_APPS))

