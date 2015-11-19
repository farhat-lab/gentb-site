"""Production settings and globals."""

from __future__ import absolute_import
import os
from os.path import join, normpath, dirname, isfile
from .base import *

# Normally you should not import ANYTHING from Django directly
# into your settings, but ImproperlyConfigured is an exception.
from django.core.exceptions import ImproperlyConfigured


#def get_env_setting(setting):
#    """ Get the environment setting or return exception """
#    try:
#        return environ[setting]
#    except KeyError:
#        error_msg = "Set the %s env variable" % setting
#        raise ImproperlyConfigured(error_msg)
import json

IS_HTTPS_SITE = True  # Used for callback urls when request object not available

########## CONFIGURATION FROM JSON FILE

json_secrets_fname = join( dirname(abspath(__file__)), "secret_settings_prod_hms.json")
if not isfile(json_secrets_fname):
    raise ValueError('JSON file in settings does not exist: %s' % json_secrets_fname)
JSON_SECRETS = json.loads(open(json_secrets_fname, 'r').read())
try:
    JSON_SECRETS = json.loads(open(json_secrets_fname, 'r').read())
except:
    raise Exception("Failed to parse JSON file for settings: %s" % json_secrets_fname)

########## END CONFIGURATION FROM JSON FILE
USE_X_FORWARDED_HOST = True

MEDIA_URL = '/tb/media/'
STATIC_URL = '/tb/static/'

STATIC_ROOT = '/www/gentb.hms.harvard.edu/docroot/tb/static'
MEDIA_ROOT = '/www/gentb.hms.harvard.edu/docroot/tb/media'
ROOT_URLCONF = '%s.urls_prod' % SITE_NAME

DEBUG = True

########## ADMIN CONTACTS
ADMINS = JSON_SECRETS['ADMINS']
TB_ADMINS = JSON_SECRETS['TB_ADMINS']
########## END ADMIN CONTACTS


########## HOST CONFIGURATION
# See: https://docs.djangoproject.com/en/1.5/releases/1.5/#allowed-hosts-required-in-production
ALLOWED_HOSTS = ['gentb.hms.harvard.edu', 'gentb-app-prod01.orchestra', 'orchestraweb.hms.harvard.edu', '134.174.150.16', 'gentb-app-shared01.orchestra']
########## END HOST CONFIGURATION

########## EMAIL CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-backend
EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'

# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-host
EMAIL_HOST = JSON_SECRETS['EMAIL_SETTINGS']['EMAIL_HOST']

# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-host-user
EMAIL_HOST_USER = JSON_SECRETS['EMAIL_SETTINGS']['EMAIL_HOST_USER']

# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-host-password
EMAIL_HOST_PASSWORD = JSON_SECRETS['EMAIL_SETTINGS']['EMAIL_HOST_PASSWORD']

# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-port
EMAIL_PORT = JSON_SECRETS['EMAIL_SETTINGS']['EMAIL_PORT']

# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-subject-prefix
EMAIL_SUBJECT_PREFIX = '[%s] ' % SITE_NAME

DEFAULT_FROM_EMAIL = JSON_SECRETS['EMAIL_SETTINGS']['DEFAULT_FROM_EMAIL']
# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-use-tls
EMAIL_USE_TLS = JSON_SECRETS['EMAIL_SETTINGS']['EMAIL_USE_TLS']

# See: https://docs.djangoproject.com/en/dev/ref/settings/#server-email
SERVER_EMAIL = EMAIL_HOST_USER
########## END EMAIL CONFIGURATION

########## DATABASE CONFIGURATION
DATABASES = {
    'default': {
        'ENGINE': JSON_SECRETS['DATABASE_SETTINGS']['ENGINE'],
        'NAME': JSON_SECRETS['DATABASE_SETTINGS']['NAME'],
        'USER': JSON_SECRETS['DATABASE_SETTINGS']['USER'],
        'PASSWORD': JSON_SECRETS['DATABASE_SETTINGS']['PASSWORD'],
        'HOST': JSON_SECRETS['DATABASE_SETTINGS']['HOST'],
        'PORT': JSON_SECRETS['DATABASE_SETTINGS']['PORT'],
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


########## SECRET CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
SECRET_KEY = JSON_SECRETS['SECRET_KEY']
#get_env_setting('SECRET_KEY')
########## END SECRET CONFIGURATION


########## TB UPLOADED DATAFILE DIRECTORY

TB_SHARED_DATAFILE_DIRECTORY = JSON_SECRETS['TB_SHARED_DATAFILE_DIRECTORY']
if not os.path.isdir(TB_SHARED_DATAFILE_DIRECTORY):
    raise Exception('Directory for uploaded TB files doesn\'t exist: %s' % TB_SHARED_DATAFILE_DIRECTORY)

########## END TB UPLOADED DATAFILE DIRECTORY

########## DROPBOX_ACCESS_TOKEN

DROPBOX_ACCESS_TOKEN = JSON_SECRETS['DROPBOX_ACCESS_TOKEN']
# see This functionality uses the Dropbox Core API to retrieve metadata from a shared link.
#    https://blogs.dropbox.com/developers/2015/08/new-api-endpoint-shared-link-metadata/

########## END DROPBOX_ACCESS_TOKEN

ERROR_LOG_FILEPATH = '/www/gentb.hms.harvard.edu/logging/gentb-django/gentb.log'
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
            'filename': ERROR_LOG_FILEPATH,
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
