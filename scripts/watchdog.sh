#!/bin/bash
#
# Watch for a refresh file and restart server.
#
if [[ -f restart_me ]]; then
  rm restart_me
  ./scripts/supervisorctl.sh restart gentb_website
fi
