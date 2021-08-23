from django.contrib.admin import ModelAdmin, StackedInline, register

from .models import CustomMap, MapRow

@register(CustomMap)
class CustomMapAdmin(ModelAdmin):
    list_display = ('name',)
    search_fields = ('name',)

@register(MapRow)
class MapRowAdmin(ModelAdmin):
    list_display = ('parent_map', 'country', 'drug')
    list_filter = ('parent_map', 'drug')
    search_fields = ('parent_map__name', 'country__name')
