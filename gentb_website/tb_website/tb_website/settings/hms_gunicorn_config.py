import os, sys
sys.path.append('/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website')
sys.path.append('/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website')
#os.chdir('/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/tb_website')

bind = '0.0.0.0:9001'
pidfile=None #'gunicorn_pid'

workers = 1
#worker_class = 'sync'
#worker_connections = 1000
timeout = 30
keepalive = 2

#logfile=/www/gentb.hms.harvard.edu/logging/gunicorn
#errorlog = '-'
loglevel = 'info'
accesslog = '/www/gentb.hms.harvard.edu/logging/gunicorn/gunicorn_access.log'
errorlog = '/www/gentb.hms.harvard.edu/logging/gunicorn/gunicorn_error.log'

proc_name = 'gunicorn_gentb' # default: None
