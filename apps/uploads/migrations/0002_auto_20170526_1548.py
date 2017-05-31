# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('uploads', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DropboxUploadFile',
            fields=[
                ('uploadfile_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='uploads.UploadFile')),
                ('url', models.URLField()),
            ],
            bases=('uploads.uploadfile',),
        ),
        migrations.AlterModelTable(
            name='uploadfile',
            table=None,
        ),
    ]
