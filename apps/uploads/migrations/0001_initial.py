# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='UploadFile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.SlugField(max_length=128)),
                ('filename', models.CharField(max_length=255, null=True)),
                ('file_directory', models.CharField(max_length=255)),
                ('size', models.PositiveIntegerField(null=True)),
                ('icon', models.URLField(null=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('retrieval_start', models.DateTimeField(null=True, blank=True)),
                ('retrieval_end', models.DateTimeField(null=True, blank=True)),
                ('retrieval_error', models.TextField(blank=True)),
            ],
            options={
                'ordering': ('-created', 'filename'),
            },
        ),
    ]
