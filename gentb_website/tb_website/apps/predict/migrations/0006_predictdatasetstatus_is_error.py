# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0005_predictdatasetstatus_human_name'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictdatasetstatus',
            name='is_error',
            field=models.BooleanField(default=False),
            preserve_default=False,
        ),
    ]
