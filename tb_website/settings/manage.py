
MANAGER_DB = None

from .local import *

print "Manager Django Backend!"

if MANAGER_DB is not None:
    for db in DATABASES:
        if db in MANAGER_DB:
            print "Updating database with MANAGER credentials."
            DATABASES[db].update(MANAGER_DB[db])

