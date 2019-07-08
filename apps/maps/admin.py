from django.contrib.admin import ModelAdmin, StackedInline, register

from .models import Place, Country, CountryHealth, CountryDetail

@register(Place)
class PlaceAdmin(ModelAdmin):
    list_display = ('name', 'country', 'rank', 'pop', 'timezone')
    search_fields = ('name',)

class HealthInline(StackedInline):
    model = CountryHealth

class DetailInline(StackedInline):
    model = CountryDetail

@register(Country)
class CountryAdmin(ModelAdmin):
    list_display = ('name', 'iso2', 'iso3')
    list_filter = ('region', 'subregion')
    search_fields = ('name', 'iso2', 'iso3', 'detail__name_short',
                     'detail__name_abbr')

    inlines = [HealthInline, DetailInline]

@register(CountryDetail)
class CountryDetailAdmin(ModelAdmin):
    list_display = ('__str__', 'name_short', 'name_abbr',
                    'rank', 'pop', 'continent')
    list_filter = ('continent', 'mapcolor')
    search_fields = ('country__name', 'name_short', 'name_abbr',)
