# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import model_utils.fields


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0001_initial'),
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
                ('dataset', models.ForeignKey(to='predict.VCFDataset')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ScriptToRun',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('name', models.CharField(max_length=100)),
                ('chosen_script', models.BooleanField(default=True)),
                ('script', models.TextField(help_text=b'Example of JSON argument: \'{"admin_url": "http://127.0.0.1:8000/tb-admin/predict/vcfdataset/3/", "callback_url": "some_url to receive results", "dataset_id": 3, "user_email": "user_who_uploaded_file@place.edu", "file1_path": ".../tb_uploaded_files/shared-files/2015/08/Predict_-_genTB_BnVjFcO.png"}\'', verbose_name=b'Command line script run by webserver.  Arguments will be passed in JSON format.')),
                ('script_args', models.TextField(help_text=b'populated on save', blank=True)),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
