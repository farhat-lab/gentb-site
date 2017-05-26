# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0008_auto_20170214_0945'),
    ]

    operations = [
        migrations.AlterField(
            model_name='predictstrain',
            name='file_one',
            field=models.ForeignKey(related_name='link_a', on_delete=django.db.models.deletion.SET_NULL, blank=True, to='uploads.UploadFile', null=True),
        ),
        migrations.AlterField(
            model_name='predictstrain',
            name='file_two',
            field=models.ForeignKey(related_name='link_b', on_delete=django.db.models.deletion.SET_NULL, blank=True, to='uploads.UploadFile', null=True),
        ),
        migrations.AlterField(
            model_name='predictstrain',
            name='piperun',
            field=models.ForeignKey(on_delete=django.db.models.deletion.SET_NULL, blank=True, to='pipeline.PipelineRun', null=True),
        ),
    ]
