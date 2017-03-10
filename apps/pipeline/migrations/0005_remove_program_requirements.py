# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pipeline', '0004_auto_20170310_1400'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='program',
            name='requirements',
        ),
    ]
