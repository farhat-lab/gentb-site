# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('uploads', '0003_resumableuploadfile'),
    ]

    operations = [
        migrations.CreateModel(
            name='ManualUploadFile',
            fields=[
                ('uploadfile_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='uploads.UploadFile', on_delete=models.CASCADE)),
                ('url', models.URLField()),
            ],
            bases=('uploads.uploadfile',),
        ),
    ]
