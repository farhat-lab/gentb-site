#!/bin/bash
#
# Watch for a refresh file and restart server.
#
if [[ -f refresh_me ]]; then
  rm refresh_me
  ./scripts/supervisorctl.sh restart all
fi
