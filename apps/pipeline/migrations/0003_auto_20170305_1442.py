# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pipeline', '0002_auto_20170217_1010'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='programrun',
            options={'ordering': ['created']},
        ),
        migrations.AddField(
            model_name='programrun',
            name='is_started',
            field=models.BooleanField(default=False),
        ),
    ]
