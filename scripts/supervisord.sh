#!/bin/bash

# Go to this directory (where the script is held)
cd "$( dirname "${BASH_SOURCE[0]}" )"
cd ..

# Set up the virtualenv
source pythonenv/bin/activate
export HOME="$PWD"

# Make double sure log directory exists
mkdir -p "$HOME/data/logs/supervisor/"

# Start supervisor
supervisord -c "$HOME/scripts/supervisord.conf" $@

