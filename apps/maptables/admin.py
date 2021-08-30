from django.contrib.admin import ModelAdmin, StackedInline, TabularInline, register

from .models import CustomMap, MapDetail, MapDataFilter, MapRow

class DetailsInline(TabularInline):
    model = MapDetail
    extra = 1

class FilterInline(StackedInline):
    model = MapDataFilter
    extra = 1

@register(CustomMap)
class CustomMapAdmin(ModelAdmin):
    list_display = ('name',)
    search_fields = ('name',)
    inlines = [DetailsInline, FilterInline]

@register(MapRow)
class MapRowAdmin(ModelAdmin):
    list_display = ('parent_map', 'country', 'drug')
    list_filter = ('parent_map', 'drug')
    search_fields = ('parent_map__name', 'country__name')
