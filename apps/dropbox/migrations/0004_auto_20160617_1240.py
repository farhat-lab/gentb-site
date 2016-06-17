# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('dropbox', '0003_auto_20160613_1306'),
    ]

    operations = [
        migrations.AddField(
            model_name='dropboxfile',
            name='filename',
            field=models.CharField(max_length=255, null=True),
        ),
        migrations.AddField(
            model_name='dropboxfile',
            name='icon',
            field=models.URLField(null=True),
        ),
        migrations.AddField(
            model_name='dropboxfile',
            name='size',
            field=models.PositiveIntegerField(null=True),
        ),
        migrations.AlterUniqueTogether(
            name='dropboxfile',
            unique_together=set([('filename', 'dataset')]),
        ),
    ]
