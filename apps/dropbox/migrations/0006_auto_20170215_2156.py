# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dropbox', '0005_auto_20170202_1444'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='dropboxfile',
            options={'ordering': ('-created', 'filename')},
        ),
        migrations.AddField(
            model_name='dropboxfile',
            name='file_directory',
            field=models.CharField(default='/tmp', max_length=255),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='dropboxfile',
            name='name',
            field=models.SlugField(max_length=128),
        ),
        migrations.AlterUniqueTogether(
            name='dropboxfile',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='dropboxfile',
            name='dataset',
        ),
        migrations.RemoveField(
            model_name='dropboxfile',
            name='fullpath',
        ),
    ]
