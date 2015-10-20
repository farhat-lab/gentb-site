# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('explore', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='exploredatafileinfo',
            options={'ordering': ('-created',), 'verbose_name': 'Explore Data File Information', 'verbose_name_plural': 'Explore Data File Information'},
        ),
    ]
