# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0003_auto_20150820_1352'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='vcfdataset',
            name='prediction_results',
        ),
    ]
