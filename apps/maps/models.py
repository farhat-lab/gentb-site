"""
Provide a place for city and country geo data for showing on maps.
"""
import os

from django.db.models import *
from django.contrib.gis.db.models import *
from django.utils.encoding import python_2_unicode_compatible

# This is where we are currently looking for data, but it could change.
url_prefix = 'http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/'

@python_2_unicode_compatible
class Country(Model):
    name = CharField(max_length=40, unique=True)
    name_short = CharField(max_length=36, null=True, blank=True)
    name_abbr = CharField(max_length=13, null=True, blank=True)

    iso2 = CharField(max_length=5, unique=True, db_index=True)
    iso3 = CharField(max_length=5, unique=True, db_index=True)

    rank = FloatField(null=True, blank=True)
    mapcolor = FloatField(null=True, blank=True)
    pop = FloatField(null=True, blank=True)
    gdp = FloatField(null=True, blank=True)

    continent = CharField(max_length=23, null=True, blank=True)
    region = CharField(max_length=23, null=True, blank=True)
    subregion = CharField(max_length=25, null=True, blank=True)
    geom = MultiPolygonField(srid=4326)
    objects = GeoManager()

    online_zip = os.path.join(url_prefix, 'ne_10m_admin_0_countries.zip')
    mapping_id = 'iso2'
    mapping = { 
        'iso2' : 'ISO_A2',
        'iso3' : 'ISO_A3',
        'rank' : 'scalerank',
        'name' : 'NAME_LONG',
        'name_short' : 'NAME',
        'name_abbr' : 'ABBREV',
        'mapcolor' : 'MAPCOLOR7',
        'pop' : 'POP_EST',
        'gdp' : 'GDP_MD_EST',
        'continent' : 'CONTINENT',
        'region' : 'REGION_UN',
        'subregion' : 'SUBREGION',
        'geom' : 'MULTIPOLYGON',
    }

    class Meta:
        verbose_name_plural = 'countries'
        ordering = ('-pop',)

    def __str__(self):
        return self.name

@python_2_unicode_compatible
class Place(Model):
    """A populated place from the world map source"""
    name = CharField(max_length=128)
    country = ForeignKey(Country, related_name='places')

    latitude = FloatField()
    longitude = FloatField()

    pop = FloatField()
    rank = IntegerField()

    elevation = FloatField()
    timezone = CharField(max_length=254)
    geom = MultiPointField(srid=4326)
    objects = GeoManager()

    online_zip = os.path.join(url_prefix, 'ne_10m_populated_places.zip')
    mapping_id = ('name', 'country')
    mapping = { 
        'pop' : 'GN_POP',
        'rank' : 'SCALERANK',
        'name' : 'NAME',
        'country' : {'iso2': 'ISO_A2'},
        'latitude' : 'LATITUDE',
        'longitude' : 'LONGITUDE',
        'elevation' : 'ELEVATION',
        'timezone' : 'TIMEZONE',
        'geom' : 'MULTIPOINT',
    }   

    class Meta:
        ordering = ('-pop',)

    def __str__(self):
        return self.name

