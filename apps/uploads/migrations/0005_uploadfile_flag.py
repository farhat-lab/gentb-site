# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('uploads', '0004_manualuploadfile'),
    ]

    operations = [
        migrations.AddField(
            model_name='uploadfile',
            name='flag',
            field=models.CharField(max_length=4, null=True, blank=True),
        ),
    ]
