from django.contrib import admin

from apps.predict.models import *

admin.site.register(PredictPipeline)

class PredictDatasetNoteAdmin(admin.ModelAdmin):
    list_display = ('dataset', 'title', 'modified', 'created',)
    search_fields = ('title', 'note')

admin.site.register(PredictDatasetNote, PredictDatasetNoteAdmin)

class PredictDatasetNoteInline(admin.StackedInline):
    model = PredictDatasetNote
    can_delete = True
    verbose_name_plural = 'Dataset Notes'
    extra = 0

class StrainInline(admin.StackedInline):
    model = PredictStrain
    can_delete = False
    extra = 0

class PredictDatasetAdmin(admin.ModelAdmin):
    inlines = (StrainInline, PredictDatasetNoteInline)
    save_on_top = True
    search_fields = ('title', 'user__first_name', 'user__last_name',)
    list_display = ('title', 'user', 'get_status', 'has_prediction', 'file_directory', 'created', 'modified')
    list_filter = ['status',]
    readonly_fields = [ 'created', 'modified', 'md5', 'file_directory',
                        'user_name', 'user_email', 'user_affiliation',]

    def get_status(self, obj):
        return obj.get_status_display()

admin.site.register(PredictDataset, PredictDatasetAdmin)


