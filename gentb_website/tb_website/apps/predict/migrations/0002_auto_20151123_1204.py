# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictdataset',
            name='fastq_type',
            field=models.CharField(blank=True, help_text=b'Only used for FastQ files', max_length=50, choices=[(b'pair-ended', b'Pair-ended'), (b'single-ended', b'Single-ended')]),
        ),
        migrations.AddField(
            model_name='predictdataset',
            name='file_type',
            field=models.CharField(default='fastq', max_length=25, choices=[(b'vcf', b'VCF'), (b'fastq', b'FastQ')]),
            preserve_default=False,
        ),
    ]
