from django.contrib.admin import ModelAdmin, StackedInline, TabularInline, register

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
    list_display = ('name', 'count_rows')

    def count_rows(self, obj):
        return obj.rows.count()

@register(MapDataRow)
class MapRowAdmin(ModelAdmin):
    list_display = ('source', 'country', 'drug')
    list_filter = ('source', 'drug')
    search_fields = ('source__name', 'country__name')
