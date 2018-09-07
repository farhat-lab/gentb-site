#
# Copyright (C) 2017 Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from django.contrib import admin

try:
    from adminsortable2.admin import SortableAdminMixin, SortableInlineAdminMixin
except ImportError:
    class SortableAdminMixin(object):
        pass
    class SortableInlineAdminMixin(object):
        pass

from .models import *

class PipelineProgramInline(SortableInlineAdminMixin, admin.TabularInline):
    model = PipelineProgram


class PipelineAdmin(admin.ModelAdmin):
    filter_horizontal = ('test_files',)
    list_display = ('name', 'description', 'errors')
    inlines = (PipelineProgramInline,)

    @staticmethod
    def errors(obj):
        """Return true if the pipeline has errors (based on tests)"""
        count = 0
        for count, test in enumerate(obj.runs.filter(run_as_test__isnull=False)):
            err = test.get_errors()
            if err:
                return err
        if count == 0:
            return 'No Tests!'
        return ''


admin.site.register(Pipeline, PipelineAdmin)

class ProgramAdmin(admin.ModelAdmin):
    list_display = ('name', 'description', 'keep',)
    filter_horizontal = ('files', 'test_files')

admin.site.register(Program, ProgramAdmin)

admin.site.register(ProgramFile)

class ProgramRunInline(admin.TabularInline):
    fields = ('program', 'job_id', 'is_submitted', 'is_started', 'is_complete', 'is_error', 'debug_text', 'error_text')
    model = ProgramRun
    extra = 0

class PipelineRunAdmin(admin.ModelAdmin):
    actions = ['all_stop']
    list_display = ('name', 'created', 'pipeline', 'status', 'age')
    search_fields = ['name', 'pipeline__name']
    list_filter = ['pipeline']
    inlines = (ProgramRunInline,)

    def status(self, obj):
        for prog in obj.programs.all():
            if prog.is_complete:
                continue
            if prog.is_error:
                return 'Error: %s' % str(prog.program)
            if prog.is_started:
                return 'Running: %s' % str(prog.program)
            if prog.is_submitted:
                return 'Waiting: %s' % str(prog.program)
        return "Complete"

    def age(self, obj):
        if obj.modified and obj.created:
            return obj.modified - obj.created
        return '-'

    def all_stop(modeladmin, request, queryset):
        for run in queryset.all():
            run.stop_all(msg='Admin Stopped this Program')
    all_stop.short_description = "Emergency All Stop"

admin.site.register(PipelineRun, PipelineRunAdmin)
