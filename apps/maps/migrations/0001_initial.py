# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from ..gis import MultiPolygonField, MultiPointField


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Country',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=40)),
                ('name_short', models.CharField(max_length=36, null=True, blank=True)),
                ('name_abbr', models.CharField(max_length=13, null=True, blank=True)),
                ('iso2', models.CharField(unique=True, max_length=5, db_index=True)),
                ('iso3', models.CharField(unique=True, max_length=5, db_index=True)),
                ('rank', models.FloatField(null=True, blank=True)),
                ('mapcolor', models.FloatField(null=True, blank=True)),
                ('pop', models.FloatField(null=True, blank=True)),
                ('gdp', models.FloatField(null=True, blank=True)),
                ('continent', models.CharField(max_length=23, null=True, blank=True)),
                ('region', models.CharField(max_length=23, null=True, blank=True)),
                ('subregion', models.CharField(max_length=25, null=True, blank=True)),
                ('geom', MultiPolygonField(srid=4326)),
            ],
        ),
        migrations.CreateModel(
            name='Place',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=128)),
                ('latitude', models.FloatField()),
                ('longitude', models.FloatField()),
                ('pop', models.FloatField()),
                ('rank', models.IntegerField()),
                ('elevation', models.FloatField()),
                ('timezone', models.CharField(max_length=254)),
                ('geom', MultiPointField(srid=4326)),
                ('country', models.ForeignKey(to='maps.Country', on_delete=models.CASCADE)),
            ],
        ),
    ]
