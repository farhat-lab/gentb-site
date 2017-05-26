# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dropbox', '0006_auto_20170215_2156'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='dropboxfile',
            name='url',
        ),
        migrations.AlterModelTable(
            name='dropboxfile',
            table='uploads_uploadedfile',
        ),
    ]
