# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0003_auto_20160525_1521'),
    ]

    operations = [
        migrations.AlterField(
            model_name='predictdatasetfile',
            name='dataset',
            field=models.ForeignKey(related_name='results', to='predict.PredictDataset', on_delete=models.CASCADE),
        ),
    ]
