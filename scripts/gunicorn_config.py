#
# Copyright 2016, see LICENSE.txt for details.
#
# pylint: disable=invalid-name
"""
Gunicorn python script, starts wsgi service with settings.
"""
import os, sys

from os.path import normpath, dirname, join

# Make sure we always know where we are when running
SCRIPT_DIR = dirname(normpath(join(os.getenv('PWD'), __file__)))
ROOT = dirname(SCRIPT_DIR)

DATA_DIR = join(ROOT, 'data')
LOG_DIR = join(DATA_DIR, 'logs', 'gunicorn')

if not os.path.isdir(LOG_DIR):
    os.makedirs(LOG_DIR)

# Add the root to the python path (so we can find mondules)
sys.path.append(ROOT)

bind = '0.0.0.0:9001'
pidfile = join(DATA_DIR, 'gunicorn.pid')

keepalive = 2
timeout = 30
workers = 1
#worker_class = 'sync'
#worker_connections = 1000

loglevel = 'info'

accesslog = join(LOG_DIR, 'access.log')
errorlog = join(LOG_DIR, 'error.log')
#errorlog = '-' # i.e. stdout

proc_name = 'gunicorn_gentb'
