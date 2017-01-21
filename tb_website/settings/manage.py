
MANAGER_DB = None

from .local import *

if MANAGER_DB is not None:
    for db in DATABASES:
        if db in MANAGER_DB:
            DATABASES[db].update(MANAGER_DB[db])

