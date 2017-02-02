# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0006_auto_20160902_1715'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='predictdatasetfile',
            name='dataset',
        ),
        migrations.DeleteModel(
            name='ScriptToRun',
        ),
        migrations.DeleteModel(
            name='PredictDatasetFile',
        ),
    ]
