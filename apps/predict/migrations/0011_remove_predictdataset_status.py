# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0010_predictdataset_delete_sources'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='predictdataset',
            name='status',
        ),
    ]
