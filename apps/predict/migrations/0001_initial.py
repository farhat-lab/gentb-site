# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import model_utils.fields


class Migration(migrations.Migration):

    dependencies = [
        ('tb_users', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DatasetScriptRun',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('notes', models.TextField(blank=True)),
                ('result_received', models.BooleanField(default=False)),
                ('result_success', models.BooleanField(default=False)),
                ('result_data', models.TextField(blank=True)),
                ('md5', models.CharField(help_text=b'auto-filled on save', max_length=40, db_index=True, blank=True)),
            ],
            options={
                'ordering': ('-modified', '-created'),
            },
        ),
        migrations.CreateModel(
            name='PipelineScriptsDirectory',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('name', models.CharField(help_text=b'helpful user name', max_length=100)),
                ('script_directory', models.TextField(help_text=b'Full path to the directory     containing the analyseNGS.pl and analyseVCF.pl pipeline scripts')),
                ('is_chosen_directory', models.BooleanField(default=True)),
            ],
            options={
                'ordering': ('is_chosen_directory', '-modified'),
                'verbose_name': 'Pipeline Scripts Directory',
                'verbose_name_plural': 'Pipeline Scripts Directory',
            },
        ),
        migrations.CreateModel(
            name='PredictDataset',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('title', models.CharField(max_length=255, verbose_name=b'Dataset title')),
                ('file_type', models.CharField(max_length=25, choices=[(b'vcf', b'VCF'), (b'fastq', b'FastQ')])),
                ('fastq_type', models.CharField(blank=True, help_text=b'Only used for FastQ files', max_length=50, choices=[(b'pair-end', b'Pair-end'), (b'single-end', b'Single-end')])),
                ('dropbox_url', models.URLField(help_text=b'https://www.dropbox.com/help/274', verbose_name=b'Dropbox link')),
                ('description', models.TextField(verbose_name=b'Dataset description')),
                ('file_directory', models.CharField(max_length=255, blank=True)),
                ('has_prediction', models.BooleanField(default=False)),
                ('md5', models.CharField(help_text=b'auto-filled on save', max_length=40, db_index=True, blank=True)),
            ],
            options={
                'ordering': ('-created', 'title'),
            },
        ),
        migrations.CreateModel(
            name='PredictDatasetFile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('name', models.CharField(help_text=b'Name of the file (w/o) the path', max_length=255)),
                ('fullpath', models.TextField(help_text=b'Full path to the file')),
                ('size', models.IntegerField(default=0, help_text=b'Size of the file in bytes')),
                ('dataset', models.ForeignKey(to='predict.PredictDataset', on_delete=models.CASCADE)),
            ],
            options={
                'ordering': ('-created', 'dataset', 'name'),
            },
        ),
        migrations.CreateModel(
            name='PredictDatasetNote',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('title', models.CharField(max_length=255)),
                ('note', models.TextField()),
                ('dataset', models.ForeignKey(to='predict.PredictDataset', on_delete=models.CASCADE)),
            ],
            options={
                'ordering': ('-modified', '-created'),
            },
        ),
        migrations.CreateModel(
            name='PredictDatasetStatus',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('human_name', models.CharField(max_length=100)),
                ('is_error', models.BooleanField()),
                ('slug', models.SlugField(blank=True)),
                ('sort_order', models.IntegerField()),
            ],
            options={
                'ordering': ('sort_order', '-name'),
                'verbose_name_plural': 'Predict Dataset Statuses',
            },
        ),
        migrations.CreateModel(
            name='ScriptToRun',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('name', models.CharField(max_length=100)),
                ('is_chosen_script', models.BooleanField(default=True)),
                ('script', models.TextField(help_text=b'Example of JSON argument: \'{"admin_url": "http://127.0.0.1:8000/tb-admin/predict/PredictDataset/3/", "callback_url": "some_url to receive results", "dataset_id": 3, "user_email": "user_who_uploaded_file@place.edu", "file1_path": ".../tb_uploaded_files/shared-files/2015/08/Predict_-_genTB_BnVjFcO.png"}\'', verbose_name=b'Command line script run by webserver.  Arguments will be passed in JSON format.')),
                ('script_args', models.TextField(help_text=b'populated on save', blank=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='predictdataset',
            name='status',
            field=models.ForeignKey(to='predict.PredictDatasetStatus', on_delete=models.CASCADE),
        ),
        migrations.AddField(
            model_name='predictdataset',
            name='user',
            field=models.ForeignKey(to='tb_users.TBUser', on_delete=models.CASCADE),
        ),
        migrations.AddField(
            model_name='datasetscriptrun',
            name='dataset',
            field=models.ForeignKey(to='predict.PredictDataset', on_delete=models.CASCADE),
        ),
    ]
