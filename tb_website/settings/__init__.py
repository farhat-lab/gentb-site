# pylint: disable=wildcard-import
"""
Find and initialize the settings/local.py - this file documents all the
settings and keys which should /NEVER/ be committed to a repository and it
seperates out the sys-admin responsibility from the programmer's.
"""

from shutil import copyfile

import logging
import sys
import os

BASE_DIR = os.path.dirname(__file__)
SETTINGS = 'local.py'

from .base import *

try:
    from .local import *
except ImportError:
    TARGET = os.path.join(BASE_DIR, SETTINGS)
    if not os.path.exists(TARGET):
        for template in (TARGET + '.template', TARGET[:-3] + '.template'):
            if os.path.exists(template):
                copyfile(template, TARGET)
                break
    try:
        from .local import *
    except ImportError:
        pass

for n in list(globals())[:]:
    v = globals()[n]
    if n.split('_')[-1] in ('DIR', 'DIRECTORY', 'PATH', 'ROOT'):
        if 'FORMAT' in n or 'LIBRARY' in n:
            continue
        if v and v[0] != '/':
            continue
        if not os.path.isdir(v):
            try:
                sys.stderr.write("INFO: Making directory: %s for %s\n" % (v, n))
                os.makedirs(v)
            except:
                if not v.startswith('/n/groups'):
                    raise IOError("Failed to make directory: %s" % v)

