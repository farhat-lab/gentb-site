#!/usr/bin/env python
#
# Companion to shell.py, script to watch for a returned value
#

import os
import sys
import time

w = open('/tmp/watch.log', 'a')

def wait(ret, cadence=0.01, *others):
    w.write("WAIT: %s (%0.2f)\n" % (ret, float(cadence)))
    while not os.path.isfile(ret):
        w.write("0\n")
        time.sleep(float(cadence))
    w.write("1\n")
    with open(ret, 'r') as fhl:
        if int(fhl.read().strip()) != 0:
            sys.exit(1)

if __name__ == '__main__':
    wait(*sys.argv[1:])
