from django.contrib.admin import *

from .models import Place, Country

class PlaceAdmin(ModelAdmin):
    list_display = ('name', 'country', 'rank', 'pop', 'timezone')
    search_fields = ('name',)

site.register(Place, PlaceAdmin)


class CountryAdmin(ModelAdmin):
    list_display = ('name', 'iso2', 'iso3', 'name_short', 'name_abbr', 'rank', 'pop', 'continent')
    list_filter = ('continent', 'region')
    search_fields = ('name',)

site.register(Country, CountryAdmin)

