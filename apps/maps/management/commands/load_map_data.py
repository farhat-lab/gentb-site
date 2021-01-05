"""
Load the country and city data
"""

import sys
import os

from zipfile import ZipFile, BadZipFile

import requests

from django.conf import settings
from django.core.management.base import BaseCommand
from django.contrib.gis.utils import LayerMapping

from apps.maps.models import Country, CountryDetail, Place
from apps.mutations.utils import StatusBar

class FilteredMapping(LayerMapping):
    def feature_kwargs(self, feat):
        def filter_value(value):
            return None if value in (-99, -99.0, '-99') else value

        def get_first_value(*keys, **kw):
            for key in keys:
                value = filter_value(feat[key].value)
                if value is not None:
                    return value
            if 'default' in kw:
                return kw['default']
            if kw.get('required', True):
                raise ValueError("No non-null value available")

        kwargs = super(FilteredMapping, self).feature_kwargs(feat)

        # Some of the values are reported as -99, which is wrong
        for key, value in kwargs.items():
            kwargs[key] = filter_value(value)

        if 'iso2' in kwargs:
            if kwargs['iso3'] is None:
                kwargs['iso3'] = get_first_value('WB_A3', 'ADM0_A3_US', 'BRK_A3')
            if kwargs['iso2'] is None:
                kwargs['iso2'] = get_first_value('WB_A2', default=kwargs['iso3'][2:])

        # The gis mapping module doesn't update the map, it adds points to it.
        # which leads to a VERY large map full of duplicate points. Clean the
        # existing data with the expectation that the geom will be repopulated.
        if self.unique:
            u_kwargs = self.unique_kwargs(kwargs)
            self.model.objects.using(self.using).filter(**u_kwargs).delete()

        return kwargs


class Command(BaseCommand):
    help = __doc__
    DATA_DIR = os.path.join(settings.DATA_ROOT, 'maps')

    def download(self, url, force=False):
        """Download the given url"""
        save_as, filename = self.url_to_filename(url)

        if os.path.isfile(save_as) and not force:
            list(StatusBar(filename, 0, []))
            return save_as

        response = requests.get(url, stream=True)
        with open(save_as, 'wb') as fhl:
            for chunk in StatusBar(filename,
              response.headers.get('content-length'),
              response.iter_content(chunk_size=4096)):
                fhl.write(chunk)

        return save_as

    def url_to_filename(self, url):
        """Where would I save the data?"""
        filename = os.path.basename(url)
        return os.path.join(self.DATA_DIR, filename), filename

    def handle(self, verbose=True, **options):
        if not os.path.isdir(self.DATA_DIR):
            os.makedirs(self.DATA_DIR)

        for model in (Country, CountryDetail, Place):
            if not hasattr(model, 'online_zip'):
                continue

            try:
                zfile = ZipFile(self.download(model.online_zip))
            except BadZipFile:
                try:
                    fileout = self.download(model.online_zip, force=True)
                    zfile = ZipFile(fileout)
                except BadZipFile:
                    sys.stderr.write(f"Could not download: {model.online_zip}\n")
                    sys.stderr.write(f"Please manually download, confirm it's a"\
                                     f"zip file and save to {fileout}")
                    break

            shp = None
            for filename in zfile.namelist():
                if filename.rsplit('.', 1)[-1] in ('prj', 'shp', 'shx', 'dbf'):
                    zfile.extract(filename, path=self.DATA_DIR)
                if filename.endswith('.shp'):
                    shp = os.path.join(self.DATA_DIR, filename)

            if shp is None or not os.path.isfile(shp):
                sys.stderr.write("Shape file %s is not found\n" % shp)
                continue

            lm = FilteredMapping(model, shp, model.mapping,
                    transform=True, unique=model.mapping_id)
            lm.save(strict=False, verbose=True)

