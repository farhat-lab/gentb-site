"""
Provide a place for city and country geo data for showing on maps.
"""
import os

from django.db.models import (
    Model, OneToOneField, ForeignKey, IntegerField, FloatField, CharField, CASCADE,
)

from .gis import GeoManager, MultiPolygonField, MultiPointField

# This is where we are currently looking for data, but it could change.
URL_PREFIX = 'http://www.naturalearthdata.com/'\
    'http//www.naturalearthdata.com/download/10m/cultural/'

# UN Regions and sub-regions
REGIONS = (
    (0, 'Other'),
    (1, 'World'),
    (2, 'Agrica'),
    (9, 'Oceania'),
    (19, 'Americas'),
    (21, 'North America'),
    (142, 'Asia'),
    (150, 'Europe'),
    (419, 'Latin America and the Caribbean'),
)
SUB_REGIONS = (
    (14, 'Eastern Africa'),
    (17, 'Middle Africa'),
    (15, 'Northern Africa'),
    (18, 'Southern Africa'),
    (11, 'Western Africa'),

    (29, 'Caribbean'),
    (13, 'Central America'),
    (5, 'South America'),

    (143, 'Central Asia'),
    (30, 'Eastern Asia'),
    (34, 'Southern Asia'),
    (35, 'South-Eastern Asia'),
    (145, 'Western Asia'),

    (151, 'Eastern Europe'),
    (154, 'Northern Europe'),
    (39, 'Southern Europe'),
    (155, 'Western Europe'),

    (53, 'Australia and New Zealand'),
    (54, 'Melanesia'),
    (57, 'Micronesia'),
    (61, 'Polynesia'),
)

class CountryManager(GeoManager):
    def get_by_natural_key(self, cid):
        return self.get(iso2=cid)

class Country(Model):
    """Provides a simple shape and a couple of useful fields to identify a country."""
    name = CharField(max_length=128, unique=True)

    iso2 = CharField(max_length=5, unique=True, db_index=True)
    iso3 = CharField(max_length=5, unique=True, db_index=True)

    region = IntegerField(choices=REGIONS, null=True, blank=True)
    subregion = IntegerField(choices=SUB_REGIONS, null=True, blank=True)

    geom = MultiPolygonField(srid=4326)
    objects = CountryManager()

    online_zip = 'http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip'
    mapping_id = 'iso2'
    mapping = {
        'iso2': 'ISO2',
        'iso3': 'ISO3',
        'geom': 'POLYGON',
        'name': 'NAME',
        'region': 'REGION',
        'subregion': 'SUBREGION',
    }

    class Meta:
        verbose_name_plural = 'countries'
        ordering = ('name',)

    def natural_key(self):
        """We want the drug code to be the key into this table"""
        return (self.iso2,)

    def __str__(self):
        return self.name

class CountryHealth(Model):
    """Extra WHO data about a country"""
    country = OneToOneField(Country, related_name='health', on_delete=CASCADE)
    total_funding = CharField(null=True, blank=True, max_length=36, \
        help_text="Total expected funding from all sources (US Dollars)")
    hiv_incidence2018 = CharField(null=True, blank=True, max_length=36,\
        help_text="Estimated incidence of TB cases who are HIV-positive")
    household = CharField(null=True, blank=True, max_length=36,\
        help_text="Estimated average household size")
    all_tb_incidence2018 = CharField(null=True, blank=True, max_length=36,\
        help_text="Estimated incidence (all forms) per 100 000 population")
    est_mdr = FloatField(null=True, blank=True,\
        help_text="Estimated % Drug Resistance for the country")

    def __str__(self):
        return "WHO TB Data for {}".format(self.country)

class CountryDetail(Model):
    """Provides a much more detailed country outline and more fields"""
    country = OneToOneField(Country, related_name='detail', on_delete=CASCADE)

    name_short = CharField(max_length=36, null=True, blank=True)
    name_abbr = CharField(max_length=13, null=True, blank=True)

    continent = CharField(max_length=23, null=True, blank=True)

    pop = FloatField(null=True, blank=True)
    gdp = FloatField(null=True, blank=True)
    rank = FloatField(null=True, blank=True)
    mapcolor = FloatField(null=True, blank=True)

    geom = MultiPolygonField(srid=4326)

    online_zip = os.path.join(URL_PREFIX, 'ne_10m_admin_0_countries.zip')
    mapping_id = 'country'
    mapping = {
        'country' : {'iso2': 'ISO_A2'},
        'rank' : 'scalerank',
        'name_short' : 'NAME',
        'name_abbr' : 'ABBREV',
        'mapcolor' : 'MAPCOLOR7',
        'gdp' : 'GDP_MD_EST',
        'continent' : 'CONTINENT',
        'geom' : 'MULTIPOLYGON',
    }

    class Meta:
        ordering = ('-pop',)

    def __str__(self):
        return self.name_short or self.name_abbr or self.country_id


class PlaceManager(GeoManager):
    def get_by_natural_key(self, name, country):
        return self.get(name=name, country=country)

class Place(Model):
    """A populated place from the world map source"""
    name = CharField(max_length=128)
    country = ForeignKey(Country, related_name='places', on_delete=CASCADE)

    latitude = FloatField()
    longitude = FloatField()

    pop = FloatField()
    rank = IntegerField()

    elevation = FloatField()
    timezone = CharField(max_length=254)
    geom = MultiPointField(srid=4326)
    objects = PlaceManager()

    online_zip = os.path.join(URL_PREFIX, 'ne_10m_populated_places.zip')
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
        unique_together = ('name', 'country',)

    def natural_key(self):
        """Return a natural key for places"""
        return (self.name, self.country.iso2)

    def __str__(self):
        return self.name
