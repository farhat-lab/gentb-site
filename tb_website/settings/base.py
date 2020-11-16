"""Common settings and globals."""

from os.path import abspath, basename, dirname, join, normpath
from sys import path

MOD_VERSION = '0.6.0'
MOD_PACKAGE = 'gentb'

DJANGO_ROOT = normpath(join(abspath(__file__), '..', '..'))

SITE_ROOT = dirname(DJANGO_ROOT)

SITE_NAME = basename(DJANGO_ROOT)
WSGI_APPLICATION = '%s.wsgi.application' % SITE_NAME

# Absolute path where all data (media, logs, etc) should go:
DATA_ROOT = join(SITE_ROOT, 'data')

# Add our project to our pythonpath, this way we don't need to type our project
# name in our dotted import paths:
path.append(DJANGO_ROOT)

########## END PATH CONFIGURATION

DEBUG_TOOLBAR_CONFIG = {
    'JQUERY_URL': '/static/js/jquery.min.js',
}

########## DEBUG CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#debug
DEBUG = False

########## GENERAL CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#time-zone
TIME_ZONE = 'America/New_York'

# See: https://docs.djangoproject.com/en/dev/ref/settings/#language-code
LANGUAGE_CODE = 'en-us'

# See: https://docs.djangoproject.com/en/dev/ref/settings/#site-id
SITE_ID = 1

# See: https://docs.djangoproject.com/en/dev/ref/settings/#use-i18n
USE_I18N = False

# See: https://docs.djangoproject.com/en/dev/ref/settings/#use-l10n
USE_L10N = True

# See: https://docs.djangoproject.com/en/dev/ref/settings/#use-tz
USE_TZ = True

FORMAT_MODULE_PATH = [
    'tb_website.formats',
]
########## END GENERAL CONFIGURATION


########## MEDIA CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#media-root
MEDIA_ROOT = normpath(join(DATA_ROOT, 'media'))

# See: https://docs.djangoproject.com/en/dev/ref/settings/#media-url
MEDIA_URL = '/media/'
########## END MEDIA CONFIGURATION

PIPELINE_ROOT = normpath(join(DATA_ROOT, 'pipeline'))
PIPELINE_BIN = normpath(join(SITE_ROOT, 'bin'))

########## STATIC FILE CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#static-root
STATIC_ROOT = normpath(join(DATA_ROOT, 'static'))

# See: https://docs.djangoproject.com/en/dev/ref/settings/#static-url
STATIC_URL = '/static/'

# Used by resumable to store upload chunks before reconstruction
UPLOAD_ROOT = normpath(join(DATA_ROOT, 'uploads'))
UPLOAD_CACHE_ROOT = normpath(join(DATA_ROOT, 'uploads', 'url_caches'))

# See: https://docs.djangoproject.com/en/dev/ref/contrib/staticfiles/#std:setting-STATICFILES_DIRS
# Should only contain static files not included in apps
STATICFILES_DIRS = []

# See: https://docs.djangoproject.com/en/dev/ref/contrib/staticfiles/#staticfiles-finders
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
)
########## END STATIC FILE CONFIGURATION

########## FIXTURE CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#std:setting-FIXTURE_DIRS
FIXTURE_DIRS = (
    normpath(join(SITE_ROOT, 'fixtures')),
)
########## END FIXTURE CONFIGURATION


########## TEMPLATE CONFIGURATION
TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            normpath(join(SITE_ROOT, 'templates')),
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]
########## END TEMPLATE CONFIGURATION


MIDDLEWARE = [
    # Default Django middleware.
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'django.middleware.security.SecurityMiddleware',
    # App based middleware
    'tb_website.middleware.AutoBreadcrumbMiddleware',
    # Support Middleware
    'apps.versioner.middleware.VersionInformation',
]

########## URL CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#root-urlconf
ROOT_URLCONF = '%s.urls' % SITE_NAME
LOGIN_URL = 'users:login'
LOGIN_REDIRECT_URL = '/'
########## END URL CONFIGURATION


########## APP CONFIGURATION
DJANGO_APPS = (
    # Default Django apps:
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',
    'django.contrib.admin',
)

# Apps from the internet (see requirements.txt)
WEBSITE_APPS = (
    'django.contrib.gis',
    'django_spaghetti',
    'adminsortable2',
    'apps.versioner',

    # Holds dataverse file ids for Two Ravens
    'apps.explore',
)

# Apps specific for this project go here.
LOCAL_APPS = (
    'tb_website',
    'apps.tb_users',
    'apps.predict',
    'apps.pipeline.app.PipelineApp',
    'apps.uploads',
    'apps.mutations',
    'apps.maps',
)

# See: https://docs.djangoproject.com/en/dev/ref/settings/#installed-apps
INSTALLED_APPS = LOCAL_APPS + DJANGO_APPS + WEBSITE_APPS
########## END APP CONFIGURATION


########## LOGGING CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#logging
# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOG_ROOT = normpath(join(DATA_ROOT, 'logs'))
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
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': join(LOG_ROOT, 'gentb.log'),
            'maxBytes': (2 ** 20) * 5, # 5 MB
            'backupCount': 5,
            'formatter': 'verbose'
        },
        'console':{
            'level':'DEBUG',
            'class':'logging.StreamHandler',
        },
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins', 'console'],
            'level': (DEBUG and 'DEBUG' or 'ERROR'),
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

SPAGHETTI_SAUCE = {
    'apps': ['mutations', 'maps', 'predict', 'pipeline', 'uploads'],
    'show_fields': False,
    'exclude': {'auth':['user']}
}

VERSION_BRANCHES = [
    ('Incoming', 'master')
]
