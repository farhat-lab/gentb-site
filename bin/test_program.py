#!/usr/bin/env python3

import sys

for name in sys.argv[1:]:
    with open(name, 'w') as fhl:
        fhl.write(f"Created file: {name}")

