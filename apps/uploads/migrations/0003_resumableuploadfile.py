# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import apps.uploads.utils
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('uploads', '0002_auto_20170526_1548'),
    ]

    operations = [
        migrations.CreateModel(
            name='ResumableUploadFile',
            fields=[
                ('uploadfile_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='uploads.UploadFile', on_delete=models.CASCADE)),
                ('upload_id', models.SlugField(default=apps.uploads.utils.get_uuid)),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL, on_delete=models.CASCADE)),
            ],
            bases=('uploads.uploadfile',),
        ),
    ]
