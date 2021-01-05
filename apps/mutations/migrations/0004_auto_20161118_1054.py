# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mutations', '0003_migrate_drug_relationships'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='genelocus',
            unique_together={('genome', 'name')},
        ),
        migrations.AlterUniqueTogether(
            name='mutation',
            unique_together={('gene_locus', 'name')},
        ),
        migrations.RemoveField(
            model_name='genelocus',
            name='drug',
        ),
        migrations.AddField(
            model_name='mutation',
            name='old_id',
            field=models.CharField(db_index=True, max_length=50, null=True, blank=True),
        ),
    ]
