# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mutations', '0004_auto_20161118_1054'),
    ]

    operations = [
        migrations.AddField(
            model_name='mutation',
            name='predictor',
            field=models.BooleanField(default=False, help_text=b'This mutation is selected to be used in predictions and will be shown to users in the manual mutation selection process.'),
        ),
    ]
