"""
Load the country and city data
"""

import requests
import sys
import os

from zipfile import ZipFile
from StringIO import StringIO
from django.conf import settings
from django.utils import six

from django.core.management.base import BaseCommand, CommandError
from django.contrib.gis.utils import LayerMapping

from apps.maps.models import Country, CountryDetail, Place

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

    def download(self, url):
        filename = os.path.basename(url)
        save_as = os.path.join(self.DATA_DIR, filename)

        sys.stderr.write("%s X%s]\r%s [" % (filename, " " * 40, filename))
        sys.stderr.flush()

        if os.path.isfile(save_as):
            sys.stderr.write('X' * 40 + "\n")
            return save_as

        pos = 0
        down = 40
        response = requests.get(url, stream=True)
        length = response.headers.get('content-length')
        with open(save_as, 'wb') as fhl:
            for chunk in response.iter_content(chunk_size=4096):
                pos += len(chunk)
                p = int((pos / float(length)) * 40)
                if 40 - p < down:
                    sys.stderr.write('-' * (down - 40 + p))
                    down = 40 - p
                fhl.write(chunk)
        sys.stderr.write("\n")
        return save_as

    def handle(self, verbose=True, **options):
        if not os.path.isdir(self.DATA_DIR):
            os.makedirs(self.DATA_DIR)

        for model in (Country, CountryDetail, Place):
            if not hasattr(model, 'online_zip'):
                continue

            zfile = ZipFile(self.download(model.online_zip))

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

