# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pipeline', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='programrun',
            name='debug_text',
            field=models.TextField(null=True, verbose_name=b'Command and Debug', blank=True),
        ),
        migrations.AlterField(
            model_name='programrun',
            name='error_text',
            field=models.TextField(null=True, verbose_name=b'Error', blank=True),
        ),
    ]
