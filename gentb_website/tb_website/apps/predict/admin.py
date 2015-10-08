from django.contrib import admin

# Register your models here.
from apps.predict.models import VCFDatasetStatus, VCFDataset, VCFDatasetNote, DatasetScriptRun, ScriptToRun


class VCFDatasetStatusAdmin(admin.ModelAdmin):
    save_on_top = True
    list_display = ['name', 'sort_order', 'slug']
    search_fields = ['name']
admin.site.register(VCFDatasetStatus, VCFDatasetStatusAdmin)



class VCFDatasetNoteAdmin(admin.ModelAdmin):
    list_display = ('dataset', 'title', 'modified', 'created',)
    search_fields = ('title', 'note')
admin.site.register(VCFDatasetNote, VCFDatasetNoteAdmin)


class VCFDatasetNoteInline(admin.StackedInline):
    model = VCFDatasetNote
    can_delete = True
    verbose_name_plural = 'Dataset Notes'
    extra = 0


class DatasetScriptRunInline(admin.StackedInline):
    model = DatasetScriptRun
    can_delete = True
    verbose_name_plural = 'Dataset Script Runs'
    extra = 0


class VCFDatasetAdmin(admin.ModelAdmin):
    inlines = (VCFDatasetNoteInline, DatasetScriptRunInline)# VCFDatasetNoteInline, )
    save_on_top = True
    search_fields = ('title', 'user__first_name', 'user__last_name',)
    list_display = ('title', 'user', 'status', 'has_prediction', 'created', 'modified')
    list_filter = ['status', 'has_prediction']
    readonly_fields = [ 'created', 'modified', 'md5',
                        'user_name', 'user_email', 'user_affiliation', 'run_script_link']

    fieldsets = [
        (None,               {'fields': ['title', ('status', 'has_prediction')]}),
        ('Description',               {'fields': ['description']}),
        ('User',               {'fields': [('user', 'user_affiliation'), ('user_name', 'user_email')]}),
        ('Files',               {'fields': ['file1', 'file2',]}),
        ('Run Script!', {'fields': ['run_script_link']}),
        ('Timestamps/md5', {'fields': [('created', 'modified'), 'md5']}),
    ]
    '''  user = models.ForeignKey(TBUser)


    has_prediction = models.BooleanField('Has prediction results?',default=False, help_text='auto-filled on save')

    prediction_results = models.TextField( blank=True)

    md5 = models.'''
admin.site.register(VCFDataset, VCFDatasetAdmin)


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

