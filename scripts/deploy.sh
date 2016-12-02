#!/bin/bash
#
# Simple bash script for updating a live instance.
#

cd "$(dirname "$0")"
cd ..

# Fetch any online revisions and check for updates
git fetch
REVS=`git log HEAD..origin/master --oneline`

if [[ $REVS == "" ]]; then
  echo "No version"
else
  echo "NEW VERSION AVAILABLE, INSTALLING:"
  echo $REVS

  # Get a new copy of available files.
  git pull

  ./pythonenv/bin/pip install -r requirements/production.txt
  ./manage migrate
  ./manage collectstatic --noinput

  touch restart_me
fi
