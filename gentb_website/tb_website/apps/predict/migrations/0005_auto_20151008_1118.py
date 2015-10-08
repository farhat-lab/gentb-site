# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0004_remove_vcfdataset_prediction_results'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='datasetscriptrun',
            options={'ordering': ('-modified', '-created')},
        ),
    ]
