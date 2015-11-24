# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import jsonfield.fields
import model_utils.fields


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0004_auto_20151124_1543'),
    ]

    operations = [
        migrations.CreateModel(
            name='DropboxRetrievalLog',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('file_metadata', jsonfield.fields.JSONField(blank=True)),
                ('file_metadata_err_msg', models.TextField(blank=True)),
                ('selected_files', jsonfield.fields.JSONField(blank=True)),
                ('retrieval_start', models.DateTimeField(null=True, blank=True)),
                ('retrieval_end', models.DateTimeField(null=True, blank=True)),
                ('retrieval_error', models.TextField(blank=True)),
                ('fastq_pair_end_extension', models.CharField(blank=True, help_text=b'For FastQ pair-end extensions. Either "_R" or "."', max_length=20, choices=[(b'_R', b'_R'), (b'.', b'.')])),
                ('files_retrieved', models.BooleanField(default=False)),
                ('md5', models.CharField(help_text=b'auto-filled on save', max_length=40, db_index=True, blank=True)),
                ('dataset', models.OneToOneField(to='predict.PredictDataset')),
            ],
            options={
                'ordering': ('-created', 'dataset'),
            },
        ),
    ]
