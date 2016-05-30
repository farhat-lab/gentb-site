from django.contrib import admin

from .models import DropboxFile

class DropboxFileAdmin(admin.ModelAdmin):
    save_on_top = True

admin.site.register(DropboxFile, DropboxFileAdmin)
