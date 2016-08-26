#
# Simple Crontab server for gentb
#

import sys
from crontab import CronTab

from datetime import datetime

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

log = sys.stdout
sys.stdout = open('/dev/null', 'w')

# Timeout of -1 means forever
for output in cron.run_scheduler(-1, **kw):
    if output:
        log.write('\n'.join([
           ("=" * 30),
           str(datetime.now()),
           ("-" * 30),
           output]) + '\n')

