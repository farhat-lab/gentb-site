from django.contrib import admin

from .models import UploadFile

class UploadFileAdmin(admin.ModelAdmin):
    save_on_top = True

admin.site.register(UploadFile, UploadFileAdmin)
