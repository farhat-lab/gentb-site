#!/bin/bash
#
# Watch for a refresh file and restart server.
#
if [[ -f refresh_me ]]; then
  ./scripts/supervisorctl.sh restart all
  rm refresh_me
fi
