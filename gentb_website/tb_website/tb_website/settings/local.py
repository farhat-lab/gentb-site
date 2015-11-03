"""Development settings and globals."""

from __future__ import absolute_import

import os
from os.path import join, normpath, isdir, isfile
import json

from .base import *

TEST_SETUP_DIR = normpath(join(dirname(dirname(DJANGO_ROOT)), 'test_setup'))
if not isdir(TEST_SETUP_DIR):
    try:
        os.makedirs(TEST_SETUP_DIR)
    except:
        print ('Failed to create directory: %s' % TEST_SETUP_DIR)

ADMINS = (
    ('Your Name', 'your_email@example.com'),

)
TB_ADMINS = (
    ('Raman Prasad', 'raman_prasad@harvard.edu'),
)
########## DEBUG CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#debug
DEBUG = True
#DEBUG = False

# See: https://docs.djangoproject.com/en/dev/ref/settings/#template-debug
TEMPLATE_DEBUG = DEBUG
########## END DEBUG CONFIGURATION
ALLOWED_HOSTS = ('127.0.0.1',)

########## EMAIL CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-backend
EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'
#EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
########## END EMAIL CONFIGURATION


########## DATABASE CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#databases
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': normpath(join(TEST_SETUP_DIR, 'db', 'tb_website.db')),
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    }
}
########## END DATABASE CONFIGURATION


########## CACHE CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#caches
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
    }
}
########## END CACHE CONFIGURATION


########## TOOLBAR CONFIGURATION
# See: http://django-debug-toolbar.readthedocs.org/en/latest/installation.html#explicit-setup
INSTALLED_APPS += (
    'debug_toolbar',
)

MIDDLEWARE_CLASSES += (
    'debug_toolbar.middleware.DebugToolbarMiddleware',
)

DEBUG_TOOLBAR_PATCH_SETTINGS = False

# http://django-debug-toolbar.readthedocs.org/en/latest/installation.html
INTERNAL_IPS = ('127.0.0.1',)
########## END TOOLBAR CONFIGURATION


########## LOGGING CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#logging
# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false'],
            'class': 'django.utils.log.AdminEmailHandler'
        },
        'file': {
            'level': 'DEBUG',
            #'class': 'logging.FileHandler',
            'class':'logging.handlers.RotatingFileHandler',
            'filename': normpath(join(TEST_SETUP_DIR, 'gentb.log')),
            'maxBytes': 1024*1024*5, # 5 MB
            'backupCount': 5,
            'formatter': 'verbose'
        },
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
        },
        'django': {
            'handlers':['file'],
            'propagate': True,
            'level':'ERROR',
        },
        'apps': {
            'handlers': ['file'],
            'propagate': True,
            'level': 'DEBUG',
        },
    },
    'formatters': {
        'verbose': {
            'format' : "[%(asctime)s] %(levelname)s [%(name)s:%(lineno)s] %(message)s",
            'datefmt' : "%d/%b/%Y %H:%M:%S"
        },
        'simple': {
            'format': '%(levelname)s %(message)s'
        },
    },
}
########## END LOGGING CONFIGURATION

########## TB UPLOADED DATAFILE DIRECTORY

TB_SHARED_DATAFILE_DIRECTORY = normpath(join(TEST_SETUP_DIR, 'tb_uploaded_files'))
if not isdir(TB_SHARED_DATAFILE_DIRECTORY):
    try:
        os.makedirs(TB_SHARED_DATAFILE_DIRECTORY)
    except:
        print ('Failed to create directory: %s' % TB_SHARED_DATAFILE_DIRECTORY)

########## END TB UPLOADED DATAFILE DIRECTORY


########## EMAIL SETTINGS

JSON_SETTINGS_SECRETS_FNAME = join(dirname(abspath(__file__)), "secret_settings_local.json")
if not isfile(JSON_SETTINGS_SECRETS_FNAME):
    raise Exception('JSON_SETTINGS_SECRETS_FNAME JSON file not found: %s' % JSON_SETTINGS_SECRETS_FNAME)

try:
    JSON_SETTINGS_SECRETS = json.loads(open(JSON_SETTINGS_SECRETS_FNAME, 'r').read())
except:
    raise Exception('Could not parse worldmap_secrets_fname JSON file: %s' % JSON_SETTINGS_SECRETS_FNAME)

EMAIL_HOST = JSON_SETTINGS_SECRETS['EMAIL_HOST']
EMAIL_PORT = JSON_SETTINGS_SECRETS['EMAIL_PORT']
EMAIL_HOST_USER = JSON_SETTINGS_SECRETS['EMAIL_HOST_USER']
EMAIL_HOST_PASSWORD = JSON_SETTINGS_SECRETS['EMAIL_HOST_PASSWORD']
DEFAULT_FROM_EMAIL = JSON_SETTINGS_SECRETS['DEFAULT_FROM_EMAIL']
EMAIL_USE_TLS = JSON_SETTINGS_SECRETS['EMAIL_USE_TLS']
########## END EMAIL SETTINGS

########## DROPBOX_ACCESS_TOKEN

DROPBOX_ACCESS_TOKEN = JSON_SETTINGS_SECRETS['DROPBOX_ACCESS_TOKEN']
# see This functionality uses the Dropbox Core API to retrieve metadata from a shared link.
#    https://blogs.dropbox.com/developers/2015/08/new-api-endpoint-shared-link-metadata/

########## END DROPBOX_ACCESS_TOKEN
