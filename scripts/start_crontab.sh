#!/bin/bash

# Make sure supervisord's stop command get's passed on python
trap "kill -- -$$" EXIT

cd "$( dirname "${BASH_SOURCE[0]}" )"
cd ..

source pythonenv/bin/activate

source /opt/lsf/conf/profile.lsf

python "$PWD/scripts/crontab_server.py" $@ 
 
