from django.contrib import admin

from apps.predict.models import *

class PredictPipelineAdmin(admin.ModelAdmin):
    list_display = ('pipeline', 'is_default', 'is_retired')
    search_fields = ('pipeline__name', 'pipeline__description')

admin.site.register(PredictPipeline, PredictPipelineAdmin)

class PredictDatasetNoteAdmin(admin.ModelAdmin):
    list_display = ('dataset', 'title', 'modified', 'created',)
    search_fields = ('title', 'note')

admin.site.register(PredictDatasetNote, PredictDatasetNoteAdmin)

class ResultsInline(admin.TabularInline):
    raw_id_fields = ('drug',)
    model = PredictResult
    can_delete = False
    extra = 0

class PredictStrainAdmin(admin.ModelAdmin):
    inlines = (ResultsInline,)
    actions = ['clear_result_cache']
    raw_id_fields = ('piperun', 'pipeline', 'dataset', 'file_one', 'file_two')
    search_fields = ("name", "dataset__title")
    list_display = ("name", "dataset", "file_one", "file_two", 'result_count')
    list_filter = ("pipeline",)

    def clear_result_cache(modeladmin, request, queryset):
        PredictResult.objects.filter(strain_id__in=queryset.values('pk')).delete()
    clear_result_cache.short_description = "Clear all Results Cache (regenerate)"

    def result_count(self, obj):
        return obj.results.count()

admin.site.register(PredictStrain, PredictStrainAdmin)

class PredictResultLocusInline(admin.TabularInline):
    raw_id_fields = ('result', 'locus',)
    model = PredictResultLocus
    can_delete = False
    extra = 0

class PredictResultAdmin(admin.ModelAdmin):
    inlines = (PredictResultLocusInline,)
    raw_id_field = ('strain', 'drug',)
    search_fields = ("strain__name", "strain__dataset__title")
    list_display = ("strain",)

admin.site.register(PredictResult, PredictResultAdmin)

class PredictDatasetNoteInline(admin.StackedInline):
    model = PredictDatasetNote
    can_delete = True
    verbose_name_plural = 'Dataset Notes'
    extra = 0

class StrainInline(admin.StackedInline):
    raw_id_fields = ('piperun', 'pipeline', 'dataset', 'file_one', 'file_two')
    model = PredictStrain
    can_delete = False
    extra = 0

class PredictDatasetAdmin(admin.ModelAdmin):
    inlines = (StrainInline, PredictDatasetNoteInline)
    actions = ['retry_processes']
    save_on_top = True
    search_fields = ('title', 'user__first_name', 'user__last_name',)
    list_display = ('title', 'user', 'get_status', 'has_prediction', 'has_lineages', 'file_directory', 'created', 'modified',)
    list_filter = []
    readonly_fields = [ 'created', 'modified', 'md5', 'file_directory',
                        'status', 'has_prediction', 'has_lineages', 'has_output_files',
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


