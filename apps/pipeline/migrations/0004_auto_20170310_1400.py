# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pipeline', '0003_auto_20170305_1442'),
    ]

    operations = [
        migrations.AddField(
            model_name='programrun',
            name='completed',
            field=models.DateTimeField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='programrun',
            name='started',
            field=models.DateTimeField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='programrun',
            name='submitted',
            field=models.DateTimeField(null=True, blank=True),
        ),
    ]
