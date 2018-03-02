"""
Provide access to gentb database without the website.
"""

import os

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "tb_website.settings.gentb")

import django
django.setup()

from apps.maps.models import *
from apps.predict.models import *
from apps.mutations.models import *
from apps.uploads.models import *
from apps.pipeline.models import *

from django.db import connections
from django.db.utils import OperationalError

db_conn = connections['default']
try:
    db_cursor = db_conn.cursor()
except OperationalError as err:
    sys.stderr.write("\nCan't connect to database: {}\n\n".format(err[1]))

from django.conf import settings
if settings.DB_IS_FILE:
    print("Making sure database schema is up to date.")
    from django.core.management import execute_from_command_line
    execute_from_command_line(['', 'migrate'])

