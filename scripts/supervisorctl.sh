#!/bin/bash

# Go to this directory (where the script is held)
cd "$( dirname "${BASH_SOURCE[0]}" )"
cd ..

# Set up the virtualenv
source pythonenv/bin/activate
export HOME="$PWD"

supervisorctl -c "$HOME/scripts/supervisord.conf" $@

