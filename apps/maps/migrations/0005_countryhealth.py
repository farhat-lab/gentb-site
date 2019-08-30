# -*- coding: utf-8 -*-
# Generated by Django 1.11.20 on 2019-07-08 08:42
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('maps', '0004_auto_20181107_1309'),
    ]

    operations = [
        migrations.CreateModel(
            name='CountryHealth',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('est_mdr', models.FloatField(blank=True, help_text=b'Estimated % Drug Resistance for the country', null=True)),
                ('country', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='health', to='maps.Country')),
            ],
        ),
    ]