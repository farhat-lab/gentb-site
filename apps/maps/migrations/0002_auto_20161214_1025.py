# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from ..gis import MultiPolygonField


class Migration(migrations.Migration):

    dependencies = [
        ('maps', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='CountryDetail',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name_short', models.CharField(max_length=36, null=True, blank=True)),
                ('name_abbr', models.CharField(max_length=13, null=True, blank=True)),
                ('continent', models.CharField(max_length=23, null=True, blank=True)),
                ('pop', models.FloatField(null=True, blank=True)),
                ('gdp', models.FloatField(null=True, blank=True)),
                ('rank', models.FloatField(null=True, blank=True)),
                ('mapcolor', models.FloatField(null=True, blank=True)),
                ('geom', MultiPolygonField(srid=4326)),
            ],
            options={
                'ordering': ('-pop',),
            },
        ),
        migrations.AlterModelOptions(
            name='country',
            options={'ordering': ('name',), 'verbose_name_plural': 'countries'},
        ),
        migrations.AlterField(
            model_name='country',
            name='name',
            field=models.CharField(unique=True, max_length=128),
        ),
        migrations.AlterModelOptions(
            name='place',
            options={'ordering': ('-pop',)},
        ),
        migrations.RemoveField(
            model_name='country',
            name='continent',
        ),
        migrations.RemoveField(
            model_name='country',
            name='gdp',
        ),
        migrations.RemoveField(
            model_name='country',
            name='mapcolor',
        ),
        migrations.RemoveField(
            model_name='country',
            name='name_abbr',
        ),
        migrations.RemoveField(
            model_name='country',
            name='name_short',
        ),
        migrations.RemoveField(
            model_name='country',
            name='pop',
        ),
        migrations.RemoveField(
            model_name='country',
            name='rank',
        ),
        migrations.RemoveField(
            model_name='country',
            name='region',
        ),
        migrations.RemoveField(
            model_name='country',
            name='subregion',
        ),
        migrations.AlterField(
            model_name='place',
            name='country',
            field=models.ForeignKey(related_name='places', to='maps.Country', on_delete=models.CASCADE),
        ),
        migrations.AddField(
            model_name='countrydetail',
            name='country',
            field=models.OneToOneField(related_name='detail', to='maps.Country', on_delete=models.CASCADE),
        ),
        migrations.AddField(
            model_name='country',
            name='region',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='country',
            name='subregion',
            field=models.IntegerField(null=True, blank=True),
        ),
    ]
