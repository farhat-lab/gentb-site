from django.contrib.admin import *

from .models import Place, Country, CountryDetail

class PlaceAdmin(ModelAdmin):
    list_display = ('name', 'country', 'rank', 'pop', 'timezone')
    search_fields = ('name',)

site.register(Place, PlaceAdmin)


class CountryAdmin(ModelAdmin):
    list_display = ('name', 'iso2', 'iso3',)
    list_filter = ('region', 'subregion')
    search_fields = ('name', 'iso2', 'iso3', 'detail__name_short', 'detail__name_abbr')

site.register(Country, CountryAdmin)


class CountryDetailAdmin(ModelAdmin):
    list_display = ('__str__', 'name_short', 'name_abbr', 'rank', 'pop', 'continent')
    list_filter = ('continent', 'mapcolor')
    search_fields = ('country__name', 'name_short', 'name_abbr',)

site.register(CountryDetail, CountryDetailAdmin)
