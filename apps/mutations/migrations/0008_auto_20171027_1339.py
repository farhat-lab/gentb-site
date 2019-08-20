# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('uploads', '0004_manualuploadfile'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('mutations', '0007_auto_20170505_1540'),
    ]

    operations = [
        migrations.AddField(
            model_name='importsource',
            name='complete',
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name='importsource',
            name='uploaded',
            field=models.ManyToManyField(to='uploads.UploadFile', blank=True),
        ),
        migrations.AddField(
            model_name='importsource',
            name='uploader',
            field=models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True, on_delete=models.CASCADE),
        ),
    ]
