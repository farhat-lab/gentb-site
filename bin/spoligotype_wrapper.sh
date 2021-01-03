#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

SPOLIGO="$DIR/spoligotype/spoligotype_info"
LOOKUP="$DIR/spoligotype/spoligotype_lookup.py"
FILENAME="$1"
SPDB="$2"

if [ ! -f "$SPOLIGO" ]; then
  cd "$DIR/spoligotype" || exit
  g++ -std=c++0x spoligotype_info.cpp -o spoligotype_info
fi

if [ "${FILENAME##*.}" == "gz" ]; then
  gunzip -c "$FILENAME" | $SPOLIGO /dev/stdin | $LOOKUP "$SPDB"
else
  $SPOLIGO "$FILENAME" | $LOOKUP "$SPDB"
fi
