from django.contrib import admin

from apps.predict.models import *
from apps.dropbox.models import DropboxFile


class PredictDatasetStatusAdmin(admin.ModelAdmin):
    save_on_top = True
    list_display = ['name', 'human_name', 'is_error', 'sort_order', 'slug']
    search_fields = ['name', 'human_name']
    list_filter = ['is_error']

admin.site.register(PredictDatasetStatus, PredictDatasetStatusAdmin)


class PredictDatasetNoteAdmin(admin.ModelAdmin):
    list_display = ('dataset', 'title', 'modified', 'created',)
    search_fields = ('title', 'note')

admin.site.register(PredictDatasetNote, PredictDatasetNoteAdmin)

class PredictDatasetFileInline(admin.StackedInline):
    model = PredictDatasetFile
    can_delete = True
    verbose_name_plural = 'Dataset Files'
    extra = 0

class PredictDatasetNoteInline(admin.StackedInline):
    model = PredictDatasetNote
    can_delete = True
    verbose_name_plural = 'Dataset Notes'
    extra = 0

class DatasetScriptRunInline(admin.StackedInline):
    model = DatasetScriptRun
    can_delete = True
    verbose_name_plural = 'Dataset Script Runs'
    extra = 0

class DropboxFileInline(admin.StackedInline):
    model = DropboxFile
    can_delete = False
    extra = 0

class PredictDatasetAdmin(admin.ModelAdmin):
    inlines = (DropboxFileInline, PredictDatasetFileInline, PredictDatasetNoteInline, DatasetScriptRunInline)
    save_on_top = True
    search_fields = ('title', 'user__first_name', 'user__last_name',)
    list_display = ('title', 'user', 'status', 'has_prediction', 'file_directory', 'created', 'modified')
    list_filter = ['status', 'has_prediction']
    readonly_fields = [ 'created', 'modified', 'md5', 'file_directory',
                        'user_name', 'user_email', 'user_affiliation',]

    fieldsets = [
        (None, {'fields': ['title', 'status', 'has_prediction']}),
        ('Description', {'fields': ['description']}),
        ('Files', {'fields': ['file_type', 'fastq_type', 'file_directory',]}),
        ('User', {'fields': [('user', 'user_affiliation'), ('user_name', 'user_email')]}),
        ('Timestamps/md5', {'fields': [('created', 'modified'), 'md5']}),
    ]

admin.site.register(PredictDataset, PredictDatasetAdmin)


class DatasetScriptRunAdmin(admin.ModelAdmin):
    save_on_top = True
    list_display = ('dataset', 'result_received', 'result_success', 'created', 'modified')
    list_filter = ['result_received', 'result_success']
    readonly_fields = [ 'created', 'modified', 'md5']

admin.site.register(DatasetScriptRun, DatasetScriptRunAdmin)

