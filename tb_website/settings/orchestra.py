"""
This module graps a configured database from 
"""

import os

from collections import defaultdict

DB_ENGINES = [
  {
    'SQLITE': 'django.db.backends.sqlite3',
    'MYSQL': 'django.db.backends.mysql',
    'PSQL': 'django.db.backends.postgresql',
    'ORACLE': None,
  },
  {
    'SQLITE': 'django.contrib.gis.db.backends.spatialite',
    'MYSQL': 'django.contrib.gis.db.backends.mysql',
    'PSQL': 'django.contrib.gis.db.backends.postgis',
    'ORACLE': None,
  },
]

def infinatedict():
    """An infinately deep defaultdict"""
    return defaultdict(infinatedict)


def parse_apache_config(filename):
    kwargs = infinatedict()
    if not os.path.isfile(filename):
        raise IOError("Apache config file is missing: %s" % filename)

    with open(filename, 'r') as fhl:
        for line in fhl.readlines():
            if line.startswith('#'):
                continue
            if line.startswith('SetEnv'):
                (name, value) = line[6:].strip().split()
                if value[0] == '"' and value[-1] == '"':
                    value = value[1:-1]
                names = name.split('_')
                target = kwargs
                for name in name.split('_'):
                    previous = target
                    target = target[name]
                previous[name] = value
    return kwargs


def get_database_config(filename, site=None, gis=False):
    config = parse_apache_config(filename)
    if not config:
        raise ValueError("Apache Config is empty: %s" % filename)

    if site is None:
        site = list(config)[0]

    if site not in config:
        raise KeyError("Apache config doesn't have site: %s" % site)

    config = config[site]
    engines = [(DB_ENGINES[gis][key], config[key]) for key, value in config.items() if key in DB_ENGINES[gis]]
    if not engines:
        raise KeyError("No supported database engine found in the apache config: %s" % filename)

    engine, db = engines[0]

    return {
        'ENGINE': engine,
        'NAME': db.get('DB', db.get('FILENAME', '')),
        'USER': db.get('USER', ''),
        'PASSWORD': db.get('PASSWORD', ''),
        'HOST': db.get('SERVER', ''),
        'PORT': db.get('PORT', ''),
    }


