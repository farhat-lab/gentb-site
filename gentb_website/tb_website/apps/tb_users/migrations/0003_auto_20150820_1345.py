# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('tb_users', '0002_tbadmincontact'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='tbadmincontact',
            options={'ordering': ('email',), 'verbose_name': 'TB Admin Contact', 'verbose_name_plural': 'TB Admin Contacts'},
        ),
    ]
