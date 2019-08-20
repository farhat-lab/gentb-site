# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Drug',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255, db_index=True)),
                ('code', models.CharField(unique=True, max_length=12, db_index=True)),
            ],
            options={
                'ordering': ('code',),
            },
        ),
        migrations.CreateModel(
            name='GeneLocus',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255, db_index=True)),
                ('drug', models.ForeignKey(related_name='gene_locuses', to='mutations.Drug', on_delete=models.CASCADE)),
            ],
            options={
                'ordering': ('name',),
            },
        ),
        migrations.CreateModel(
            name='Mutation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255, db_index=True)),
                ('order', models.IntegerField(default=0)),
                ('gene_locus', models.ForeignKey(related_name='mutations', to='mutations.GeneLocus', on_delete=models.CASCADE)),
            ],
            options={
                'ordering': ('order',),
            },
        ),
        migrations.AlterUniqueTogether(
            name='mutation',
            unique_together=set([('gene_locus', 'name')]),
        ),
        migrations.AlterUniqueTogether(
            name='genelocus',
            unique_together=set([('drug', 'name')]),
        ),
    ]
