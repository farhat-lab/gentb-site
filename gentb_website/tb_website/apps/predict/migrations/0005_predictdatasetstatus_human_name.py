# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0004_auto_20151124_1543'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictdatasetstatus',
            name='human_name',
            field=models.CharField(default='hello', max_length=100),
            preserve_default=False,
        ),
    ]
