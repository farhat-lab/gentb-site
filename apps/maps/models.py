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

class druginfo(Model):
    #def __init__(Name):
       # name = Name
        
    """can keep track of each drug info from Avika's data in each country"""
    country = OneToOneField(Country, related_name='drug', on_delete=CASCADE)
    name = CharField(null=True, blank=True,max_length=36)
    resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to {}".format(name))
    rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to {}".format(name))
    rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to {}".format(name))
    mean_resistance = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to {}".format(name))
    mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to {}".format(name))


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
    pop_dens = CharField(null=True, blank=True, max_length=36,\
        help_text="Population density (people per sq. km of land area")
    est_mdr = FloatField(null=True, blank=True,\
        help_text="Estimated % Drug Resistance for the country")
    world_bank_gdp = CharField(null=True, blank=True, max_length=36,\
        help_text="GDP (current US$)")
    #total_wealth = CharField(null=True, blank=True, max_length=36,\
    #    help_text="Total wealth per capita (constant 2014 US$)")
    #percent_resistance = CharField(null=True, blank=True,max_length=36,\
    #    help_text="The percent of isolates resistant to isonized variants")
    #confidence_interval = CharField(null=True, blank=True,max_length=36,\
    #    help_text="The confidence interval in which the percent resistance is calculated")
    sample_size = CharField(null=True, blank=True,max_length=36,\
        help_text="the number of samples from each country")
    resistant_isolates = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to rif")
    suseptable_isolates = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates suseptable to rif")
    ### AMK data
    amk_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to amk")
    amk_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to amk")
    amk_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to amk")
    amk_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to amk")
    amk_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to amk")
    #RIF data
    rif_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to rif")
    rif_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to rif")
    rif_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to rif")
    rif_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to rif")
    rif_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to rif")
    #CAP data
    cap_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to cap")
    cap_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to cap")
    cap_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to cap")
    cap_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to cap")
    cap_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to cap")
    #CIP data
    cip_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to cip")
    cip_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to cip")
    cip_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to cip")
    cip_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to cip")
    cip_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to cip")
    #EMB data
    emb_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to emb")
    emb_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to emb")
    emb_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to emb")
    emb_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to emb")
    emb_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to emb")
    #ETH data
    eth_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to eth")
    eth_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to eth")
    eth_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to eth")
    eth_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to eth")
    eth_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to eth")
    #INH data
    inh_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to inh")
    inh_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to inh")
    inh_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to inh")
    inh_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to inh")
    inh_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to inh")
    #KAN data
    kan_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to kan")
    kan_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to kan")
    kan_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to kan")
    kan_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to kan")
    kan_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to kan")
    #LEVO data
    levo_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to levo")
    levo_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to levo")
    levo_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to levo")
    levo_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to levo")
    levo_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to levo")
    #OFLX data
    oflx_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to oflx")
    oflx_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to oflx")
    oflx_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to oflx")
    oflx_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to oflx")
    oflx_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to oflx")
    #PAS data
    pas_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to pas")
    pas_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to pas")
    pas_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to pas")
    pas_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to pas")
    pas_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to pas")
    #PZA data
    pza_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to pza")
    pza_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to pza")
    pza_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to pza")
    pza_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to pza")
    pza_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to pza")
    #STR data
    str_resistant = CharField(null=True, blank=True,max_length=36,\
        help_text="Number of isolates resistant to str")
    str_rr = CharField(null=True, blank=True,max_length=36,\
        help_text="of the rif resistant isolates, isolates also resistant to str")
    str_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Of the rif suseptable isolates, isolates resistant to str")
    str_mean_rr = CharField(null=True, blank=True,max_length=36,\
         help_text="Mean resistance to str")
    str_mean_rs = CharField(null=True, blank=True,max_length=36,\
        help_text="Mean rif suseptibility to str")


    

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
