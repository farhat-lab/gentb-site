#
# Copyright (C) 2016   Dr. Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the 
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
A basic set of administrative functions for genes and mutations
"""

from django.contrib.admin import *
from .models import *
from .forms import DrugForm

class ImportSourceAdmin(ModelAdmin):
    list_display = ('name', 'created', 'uploader', 'complete')
    readonly_fields = ('uploaded',)

site.register(ImportSource, ImportSourceAdmin)

class DrugAdmin(ModelAdmin):
    list_display = ('__str__', 'mutation_count')

    def get_form(self, *args, **kw):
        return DrugForm

    def mutation_count(self, obj):
        return obj.mutations.count()

site.register(DrugClass)
site.register(Drug, DrugAdmin)

class MutationInline(StackedInline):
    model = Mutation
    extra = 2

class GeneLocusAdmin(ModelAdmin):
    inlines = (MutationInline,)
    list_filter = ('genome', 'gene_type', 'strand')
    list_display = ('name', 'description', 'previous_id', 'gene_symbol', 'mutation_count', 'genome')

    def mutation_count(self, obj):
        return obj.mutations.count()

site.register(Genome)
site.register(GeneLocus, GeneLocusAdmin)

class MutationAdmin(ModelAdmin):
    list_display = ('name', 'old_id', 'gene_locus', 'drugs_list')
    list_filter = ('predictor', 'gene_locus', 'drugs')
    search_fields = ('name', 'old_id')

    def drugs_list(self, obj):
        return ", ".join(obj.drugs.values_list('code', flat=True))

site.register(Mutation, MutationAdmin)

site.register(TargetSet)
site.register(TargetRegion)

class StrainSourceAdmin(ModelAdmin):
    list_display = ('__str__', 'old_id', 'patient_id', 'country', 'date')
    list_filter = ('importer', 'source_lab', 'patient_sex', 'patient_hiv', 'resistance_group')

site.register(StrainSource, StrainSourceAdmin)
site.register(StrainMutation)

class StrainResistanceAdmin(ModelAdmin):
    list_filter = ('resistance', 'drug')

site.register(StrainResistance, StrainResistanceAdmin)

