from django.contrib.admin import ModelAdmin, StackedInline, TabularInline, register
from django.utils.html import format_html

from .models import MapDisplayDetail, MapDisplayFilter, MapDisplay, MapDataSource, MapDataRow

class DetailInline(TabularInline):
    model = MapDisplayDetail
    extra = 1

class FilterInline(StackedInline):
    model = MapDisplayFilter
    extra = 1

@register(MapDisplay)
class MapDisplayAdmin(ModelAdmin):
    list_display = ('name', 'data', 'description')
    list_filter = ('data',)
    search_fields = ('name', 'description')
    inlines = [DetailInline, FilterInline]

@register(MapDataSource)
class MapSourceAdmin(ModelAdmin):
    list_display = ('name', 'count_rows', 'get_columns')

    def count_rows(self, obj):
        """Rows"""
        return obj.rows.count()
    count_rows.short_description = "Number of Rows"

    def get_columns(self, obj):
        return format_html("<br/>".join(obj.get_columns()))
    get_columns.short_description = "Column Names"

@register(MapDataRow)
class MapRowAdmin(ModelAdmin):
    list_display = ('source', 'country', 'drug')
    list_filter = ('source', 'drug')
    search_fields = ('source__name', 'country__name')
