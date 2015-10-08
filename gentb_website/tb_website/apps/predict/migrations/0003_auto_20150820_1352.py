# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0002_datasetscriptrun_scripttorun'),
    ]

    operations = [
        migrations.RenameField(
            model_name='scripttorun',
            old_name='chosen_script',
            new_name='is_chosen_script',
        ),
    ]
