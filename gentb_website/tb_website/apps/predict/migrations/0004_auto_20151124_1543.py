# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0003_pipelinescriptsdirectory'),
    ]

    operations = [
        migrations.AlterField(
            model_name='predictdataset',
            name='fastq_type',
            field=models.CharField(blank=True, help_text=b'Only used for FastQ files', max_length=50, choices=[(b'pair-end', b'Pair-end'), (b'single-end', b'Single-end')]),
        ),
    ]
