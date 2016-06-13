# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('dropbox', '0002_auto_20160613_1047'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='dropboxfile',
            unique_together=set([('name', 'dataset')]),
        ),
    ]
