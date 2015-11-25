from django.contrib import admin

# Register your models here.
from apps.predict.models import PredictDatasetStatus, PredictDataset,\
            PredictDatasetNote, DatasetScriptRun,\
            ScriptToRun, PipelineScriptsDirectory
from apps.dropbox_helper.models import DropboxRetrievalLog

from apps.predict.admin_forms import PipelineScriptsDirectoryForm


class PredictDatasetStatusAdmin(admin.ModelAdmin):
    save_on_top = True
    list_display = ['name', 'human_name', 'sort_order', 'slug']
    search_fields = ['name', 'human_name']
admin.site.register(PredictDatasetStatus, PredictDatasetStatusAdmin)


class PredictDatasetNoteAdmin(admin.ModelAdmin):
    list_display = ('dataset', 'title', 'modified', 'created',)
    search_fields = ('title', 'note')
admin.site.register(PredictDatasetNote, PredictDatasetNoteAdmin)


class PipelineScriptsDirectoryAdmin(admin.ModelAdmin):
    form = PipelineScriptsDirectoryForm
    save_on_top = True
    list_display = ['name', 'is_chosen_directory', 'modified', 'script_directory']
    search_fields = ['name']
admin.site.register(PipelineScriptsDirectory, PipelineScriptsDirectoryAdmin)


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


"""
class DropboxDataSourceAdmin(admin.ModelAdmin):
    save_on_top = True
    list_display = ('dataset', 'files_retrieved', 'created', 'dropbox_url')
    list_filter = ('files_retrieved', )
admin.site.register(DropboxDataSource, DropboxDataSourceAdmin)

class DropboxDataSourceInline(admin.StackedInline):
    model = DropboxDataSource
    can_delete = True
    extra = 0
"""

class DropboxDataSourceInline(admin.StackedInline):
    model = DropboxRetrievalLog
    can_delete = False
    extra = 0

class PredictDatasetAdmin(admin.ModelAdmin):
    inlines = (DropboxDataSourceInline, PredictDatasetNoteInline, DatasetScriptRunInline)# PredictDatasetNoteInline, )
    save_on_top = True
    search_fields = ('title', 'user__first_name', 'user__last_name',)
    list_display = ('title', 'user', 'status', 'has_prediction', 'created', 'modified')
    list_filter = ['status', 'has_prediction']
    readonly_fields = [ 'created', 'modified', 'md5', 'file_directory',
                        'user_name', 'user_email', 'user_affiliation',]

    fieldsets = [
        (None,               {'fields': ['title', 'status', 'has_prediction']}),
        ('Description',               {'fields': ['description']}),
        ('Files',               {'fields': ['file_type', 'fastq_type', 'dropbox_url', 'file_directory',]}),
        ('User',               {'fields': [('user', 'user_affiliation'), ('user_name', 'user_email')]}),
        #('Run Script!', {'fields': ['run_script_link']}),
        ('Timestamps/md5', {'fields': [('created', 'modified'), 'md5']}),
    ]

admin.site.register(PredictDataset, PredictDatasetAdmin)


class DatasetScriptRunAdmin(admin.ModelAdmin):
    save_on_top = True
    list_display = ('dataset', 'result_received', 'result_success', 'created', 'modified')
    list_filter = ['result_received', 'result_success']
    readonly_fields = [ 'created', 'modified', 'md5']
admin.site.register(DatasetScriptRun, DatasetScriptRunAdmin)

class ScriptToRunAdmin(admin.ModelAdmin):
    save_on_top = True
    list_display = ('name', 'is_chosen_script', 'script', 'modified')
    readonly_fields = ['created', 'modified', 'script_args']
    list_filter = ['is_chosen_script']
admin.site.register(ScriptToRun, ScriptToRunAdmin)
