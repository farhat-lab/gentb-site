# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mutations', '0008_auto_20171027_1339'),
    ]

    operations = [
        migrations.CreateModel(
            name='Paper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=128)),
                ('doi', models.CharField(unique=True, max_length=255)),
                ('url', models.URLField()),
                ('notes', models.TextField(null=True, blank=True)),
            ],
        ),
        migrations.AlterField(
            model_name='mutation',
            name='ecoli_aapos',
            field=models.IntegerField(help_text=b'E.coli aminoacid position used in lookups when mutations are in the genes that are in both bacteria', null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='strainsource',
            name='wgs_group',
            field=models.CharField(max_length=10, null=True, verbose_name=b'Whole Genome Sequence Group', blank=True),
        ),
        migrations.AddField(
            model_name='strainsource',
            name='source_paper',
            field=models.ForeignKey(related_name='strains', blank=True, to='mutations.Paper', null=True, on_delete=models.CASCADE),
        ),
    ]
