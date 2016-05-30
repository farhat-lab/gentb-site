# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0002_auto_20160524_0947'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='predictdataset',
            name='dropbox_url',
        ),
        migrations.AlterField(
            model_name='predictdataset',
            name='file_type',
            field=models.CharField(max_length=25, choices=[(b'vcf', b'Variant Call Format (VCF)'), (b'fastq', b'FastQ Nucleotide Sequence'), (b'manual', b'Mutations Manual Entry')]),
        ),
        migrations.AlterField(
            model_name='predictdataset',
            name='title',
            field=models.CharField(max_length=255, verbose_name=b'Dataset Title'),
        ),
    ]
