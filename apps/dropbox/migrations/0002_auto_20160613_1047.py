# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0004_auto_20160613_1047'),
        ('dropbox', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DropboxFile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('url', models.URLField()),
                ('name', models.SlugField(max_length=32)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('retrieval_start', models.DateTimeField(null=True, blank=True)),
                ('retrieval_end', models.DateTimeField(null=True, blank=True)),
                ('retrieval_error', models.TextField(blank=True)),
                ('dataset', models.ForeignKey(related_name='files', to='predict.PredictDataset')),
                ('result', models.ForeignKey(blank=True, to='predict.PredictDatasetFile', null=True)),
            ],
            options={
                'ordering': ('-created', 'dataset'),
            },
        ),
        migrations.RemoveField(
            model_name='dropboxretrievallog',
            name='dataset',
        ),
        migrations.DeleteModel(
            name='DropboxRetrievalLog',
        ),
    ]
