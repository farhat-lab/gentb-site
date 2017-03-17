# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mutations', '0005_mutation_predictor'),
    ]

    operations = [
        migrations.AddField(
            model_name='strainsource',
            name='wgs_group',
            field=models.CharField(max_length=10, null=True, verbose_name=b'Whole Gnome Sequence Group', blank=True),
        ),
    ]
