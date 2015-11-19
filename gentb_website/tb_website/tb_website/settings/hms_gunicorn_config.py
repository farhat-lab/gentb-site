
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
accesslog = '/www/gentb.hms.harvard.edu/logging/gunicorn_access.log'
errorlog = '/www/gentb.hms.harvard.edu/logging/gunicorn_error.log'

proc_name = 'gunicorn_gentb' # default: None
