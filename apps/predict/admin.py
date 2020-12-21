from django.contrib import admin

from apps.predict.models import *

admin.site.register(PredictPipeline)

class PredictDatasetNoteAdmin(admin.ModelAdmin):
    list_display = ('dataset', 'title', 'modified', 'created',)
    search_fields = ('title', 'note')

admin.site.register(PredictDatasetNote, PredictDatasetNoteAdmin)

admin.site.register(PredictStrain)

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
    #inlines = (StrainInline, PredictDatasetNoteInline)
    actions = ['retry_processes']
    save_on_top = True
    search_fields = ('title', 'user__first_name', 'user__last_name',)
    list_display = ('title', 'user', 'get_status', 'has_prediction', 'file_directory', 'created', 'modified')
    list_filter = []
    readonly_fields = [ 'created', 'modified', 'md5', 'file_directory',
                        'user_name', 'user_email', 'user_affiliation',]

    def retry_processes(modeladmin, request, queryset):
        for predict in queryset.all():
            for strain in predict.strains.all():
                if not strain.piperun:
                    continue
                for run in strain.piperun.all_programs():
                    if run.is_submitted and run.is_error and run.has_input:
                        run.is_submitted = False
                        run.error_text = ''
                        run.save()
    retry_processes.short_description = "Atempt to Retry Jobs"


admin.site.register(PredictDataset, PredictDatasetAdmin)


