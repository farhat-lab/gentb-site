# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dropbox', '0004_auto_20160617_1240'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='dropboxfile',
            name='result',
        ),
        migrations.AddField(
            model_name='dropboxfile',
            name='fullpath',
            field=models.TextField(null=True, blank=True),
        ),
    ]
