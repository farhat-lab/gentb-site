from django.contrib import admin
from apps.shared_data.models import SharedFileInfo

class SharedFileInfoAdmin(admin.ModelAdmin):
    save_on_top = True
    search_fields = ('title', 'affiliation', 'last_name', 'first_name')
    list_display = ('title', 'has_prediction', 'affiliation', 'last_name', 'first_name', 'contact_email')    
    list_filter = [ 'affiliation', 'has_prediction']
    readonly_fields = [ 'created', 'modified', 'md5', 'has_prediction']
admin.site.register(SharedFileInfo, SharedFileInfoAdmin)

