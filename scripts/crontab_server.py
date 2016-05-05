#
# Simple Crontab server for gentb
#

import sys
from crontab import CronTab

fn = 'scripts/crontab.tab'

try:
    cron = CronTab(tabfile=fn)
except IOError:
    sys.stderr.write("Can't open tabfile: %s\n" % fn)
    sys.exit(2)

kw = {
  'warp': bool('--test' in sys.argv),
}
    
if kw['warp']:
    kw['cadence'] = 1
    print "Testing mode, 1 minute == 1 second."

# Timeout of -1 means forever
for output in cron.run_scheduler(-1, **kw):
    print output

