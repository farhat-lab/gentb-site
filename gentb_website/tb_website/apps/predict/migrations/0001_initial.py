# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import django.core.files.storage
import model_utils.fields

from django.conf import settings

class Migration(migrations.Migration):

    dependencies = [
        ('tb_users', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='VCFDataset',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('title', models.CharField(max_length=255, verbose_name=b'Dataset title')),
                ('description', models.TextField(verbose_name=b'Dataset description')),
                ('file1', models.FileField(upload_to=b'shared-files/%Y/%m', storage=django.core.files.storage.FileSystemStorage(location=settings.TB_SHARED_DATAFILE_DIRECTORY), verbose_name=b'File 1')),
                ('file2', models.FileField(storage=django.core.files.storage.FileSystemStorage(location=settings.TB_SHARED_DATAFILE_DIRECTORY), upload_to=b'shared-files/%Y/%m', null=True, verbose_name=b'File 2 (optional)', blank=True)),
                ('has_prediction', models.BooleanField(default=False, help_text=b'auto-filled on save', verbose_name=b'Has prediction results?')),
                ('prediction_results', models.TextField(blank=True)),
                ('md5', models.CharField(help_text=b'auto-filled on save', max_length=40, db_index=True, blank=True)),
            ],
            options={
                'ordering': ('-created', 'title'),
                'verbose_name': 'VCF Dataset',
                'verbose_name_plural': 'VCF Datasets',
            },
        ),
        migrations.CreateModel(
            name='VCFDatasetNote',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('title', models.CharField(max_length=255)),
                ('note', models.TextField()),
                ('dataset', models.ForeignKey(to='predict.VCFDataset')),
            ],
            options={
                'ordering': ('-modified', '-created'),
            },
        ),
        migrations.CreateModel(
            name='VCFDatasetStatus',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('slug', models.SlugField(blank=True)),
                ('sort_order', models.IntegerField()),
            ],
            options={
                'ordering': ('sort_order', '-name'),
                'verbose_name': 'VCF Dataset Status',
                'verbose_name_plural': 'VCF Dataset Statuses',
            },
        ),
        migrations.AddField(
            model_name='vcfdataset',
            name='status',
            field=models.ForeignKey(to='predict.VCFDatasetStatus'),
        ),
        migrations.AddField(
            model_name='vcfdataset',
            name='user',
            field=models.ForeignKey(to='tb_users.TBUser'),
        ),
    ]
