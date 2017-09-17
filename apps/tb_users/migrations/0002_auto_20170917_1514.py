# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('tb_users', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='tbuser',
            options={'permissions': (('can_raw_upload', 'User can upload by server url.'),)},
        ),
    ]
