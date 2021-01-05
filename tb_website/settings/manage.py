
MANAGER_DB = None

from .local import *

#sys.stderr.write("Manager Django Backend!\n")

if MANAGER_DB is not None:
    for db in DATABASES:
        if db in MANAGER_DB:
            #sys.stderr.write("Updating database with MANAGER credentials.\n")
            DATABASES[db].update(MANAGER_DB[db])

