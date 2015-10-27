# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


def load_explore_links(apps, schema_editor):

    ExploreDataFileInfo = apps.get_model("explore", "ExploreDataFileInfo")

    data_file_info = ExploreDataFileInfo(name="initial entry",
            codebook_file_url="https://dataverse.harvard.edu/api/access/datafile/2694344",
            two_ravens_url="https://rserve.dataverse.harvard.edu/dataexplore/gui.html?dfId=2693726&key=c54f07b7-5098-461c-adf3-a976c0d62f6e",
            active=True)
    data_file_info.save()

class Migration(migrations.Migration):

    dependencies = [
        ('explore', '0002_auto_20151014_1324'),
    ]

    operations = [
         migrations.RunPython(load_explore_links),
    ]
