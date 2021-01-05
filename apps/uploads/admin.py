from django.contrib import admin
from django.utils.safestring import mark_safe

from .models import UploadFile

class UploadFileAdmin(admin.ModelAdmin):
    list_display = ('title', 'get_type', 'flag', 'size', 'created', 'retrieval_start', 'retrieval_end', 'retrieval_error')
    save_on_top = True

    @staticmethod
    def title(obj):
        """Returns the uploads name plus icon (if available)"""
        if obj.icon:
            return mark_safe("<img src='{0.icon}' style='width: 2em; vertical-align: middle;'> {0.name}".format(obj))
        return obj.name

admin.site.register(UploadFile, UploadFileAdmin)
