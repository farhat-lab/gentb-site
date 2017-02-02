# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0005_auto_20160617_1231'),
    ]

    operations = [
        migrations.AddField(
            model_name='datasetscriptrun',
            name='process_end',
            field=models.DateTimeField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='datasetscriptrun',
            name='process_start',
            field=models.DateTimeField(auto_now_add=True, null=True),
        ),
        migrations.AddField(
            model_name='predictdataset',
            name='status',
            field=models.PositiveIntegerField(default=1, choices=[(0, 'Dataset Deleted'), (1, 'Dataset Not Ready'), (2, 'Dataset Confirmed'), (3, 'File Retrieval Started'), (4, 'File Retrieval Failed'), (5, 'File Retrieval Success'), (6, 'Processing Started'), (7, 'Processing Success'), (8, 'Processing Failed')]),
        ),
        migrations.DeleteModel(
            name='PredictDatasetStatus',
        ),
    ]
