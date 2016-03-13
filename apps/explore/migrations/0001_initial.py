# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import model_utils.fields


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='ExploreDataFileInfo',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('name', models.CharField(help_text=b'For internal use', max_length=255)),
                ('active', models.BooleanField(default=True, help_text=b'The *most recently created* active entry will be used')),
                ('codebook_file_url', models.URLField(help_text=b'Example: https://dataverse.harvard.edu/api/access/datafile/2694344')),
                ('two_ravens_url', models.URLField(help_text=b'Example: https://rserve.dataverse.harvard.edu/dataexplore/gui.html?dfId=2693726&key=c54f07b7-5098-461c-adf3-a976c0d62f6e')),
            ],
            options={
                'ordering': ('-created',),
                'verbose_name': 'Explore Data File Information',
                'verbose_name_plural': 'Explore Data File Information',
            },
        ),
    ]
