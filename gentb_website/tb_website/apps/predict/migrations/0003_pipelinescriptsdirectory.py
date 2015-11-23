# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import model_utils.fields


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0002_auto_20151123_1204'),
    ]

    operations = [
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
    ]
