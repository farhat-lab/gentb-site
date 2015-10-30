from django.contrib import admin

from apps.dropbox_helper.models import DropboxRetrievalLog

class DropboxRetrievalLogAdmin(admin.ModelAdmin):
    save_on_top = True
admin.site.register(DropboxRetrievalLog, DropboxRetrievalLogAdmin)
