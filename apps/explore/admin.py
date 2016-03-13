from django.contrib import admin
from apps.explore.models import ExploreDataFileInfo

class ExploreDataFileInfoAdmin(admin.ModelAdmin):
    save_on_top = True
    readonly_fields = ['created', 'modified']
    list_display = ['name', 'active', 'codebook_file_url', 'two_ravens_url', 'created']
    search_fields = ['name']
    list_filter = ['active']
admin.site.register(ExploreDataFileInfo, ExploreDataFileInfoAdmin)
