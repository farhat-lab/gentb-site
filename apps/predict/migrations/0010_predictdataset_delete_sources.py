# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0009_auto_20170305_1442'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictdataset',
            name='delete_sources',
            field=models.BooleanField(default=False, help_text=b'If this is checked, we will delete all your input files downloaded from dropbox after running the predict.'),
        ),
    ]
