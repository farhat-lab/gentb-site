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

from adminsortable2.admin import SortableAdminMixin, SortableInlineAdminMixin

from .models import *

class PipelineProgramInline(SortableInlineAdminMixin, admin.TabularInline):
    model = PipelineProgram


class PipelineAdmin(admin.ModelAdmin):
    filter_horizontal = ('test_files',)
    list_display = ('name', 'description')
    inlines = (PipelineProgramInline,)


admin.site.register(Pipeline, PipelineAdmin)

class ProgramAdmin(admin.ModelAdmin):
    list_display = ('name', 'description', 'keep')
    filter_horizontal = ('files', 'test_files')

admin.site.register(Program, ProgramAdmin)

admin.site.register(ProgramFile)

class ProgramRunInline(admin.TabularInline):
    fields = ('program', 'job_id', 'is_submitted', 'is_started', 'is_complete', 'is_error', 'debug_text', 'error_text')
    model = ProgramRun
    extra = 0

class PipelineRunAdmin(admin.ModelAdmin):
    actions = ['all_stop']
    list_display = ('name', 'pipeline',)
    inlines = (ProgramRunInline,)

    def all_stop(modeladmin, request, queryset):
        for run in queryset.all():
            run.stop_all(msg='Admin Stopped this Program')
    all_stop.short_description = "Emergency All Stop"

admin.site.register(PipelineRun, PipelineRunAdmin)
