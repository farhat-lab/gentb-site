# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('maps', '0002_auto_20161214_1025'),
    ]

    operations = [
        migrations.AlterField(
            model_name='country',
            name='region',
            field=models.IntegerField(blank=True, null=True, choices=[(0, b'Other'), (1, b'World'), (2, b'Agrica'), (9, b'Oceania'), (19, b'Americas'), (21, b'North America'), (142, b'Asia'), (150, b'Europe'), (419, b'Latin America and the Caribbean')]),
        ),
        migrations.AlterField(
            model_name='country',
            name='subregion',
            field=models.IntegerField(blank=True, null=True, choices=[(14, b'Eastern Africa'), (17, b'Middle Africa'), (15, b'Northern Africa'), (18, b'Southern Africa'), (11, b'Western Africa'), (29, b'Caribbean'), (13, b'Central America'), (5, b'South America'), (143, b'Central Asia'), (30, b'Eastern Asia'), (34, b'Southern Asia'), (35, b'South-Eastern Asia'), (145, b'Western Asia'), (151, b'Eastern Europe'), (154, b'Northern Europe'), (39, b'Southern Europe'), (155, b'Western Europe'), (53, b'Australia and New Zealand'), (54, b'Melanesia'), (57, b'Micronesia'), (61, b'Polynesia')]),
        ),
    ]
