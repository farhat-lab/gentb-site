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
# pylint: disable=too-few-public-methods, missing-docstring
#
"""
A basic set of administrative functions for genes and mutations
"""

from django.contrib.admin import register, site, ModelAdmin, StackedInline, TabularInline
from .models import (
    Drug, DrugClass, DrugRegimen, ImportSource, Mutation, Genome, GeneLocus,
    GeneDrugInteraction, TargetSet, TargetRegion, StrainResistance, Paper,
    BioProject, StrainSource, StrainMutation, Lineage
)

@register(Lineage)
class LineageAdmin(ModelAdmin):
    list_display = ('slug', 'name', 'parent')

@register(ImportSource)
class ImportSourceAdmin(ModelAdmin):
    list_display = ('name', 'created', 'uploader', 'complete')
    exclude = ('uploaded',)

class InteractionInline(TabularInline):
    model = GeneDrugInteraction
    raw_id_fields = ('gene',)
    fields = ('gene', 'interaction', 'weight', 'paper')

@register(Drug)
class DrugAdmin(ModelAdmin):
    list_display = ('__str__', 'abbr', 'gene_count', 'mutation_count', 'kind', 'regimen')
    list_filter = ('kind', 'regimen')
    readonly_fields = ('mutations',)
    inlines = (InteractionInline,)
    filter_horizontal = ('loci',)

    @staticmethod
    def gene_count(obj):
        return obj.gene_interactions.count()

    @staticmethod
    def mutation_count(obj):
        return obj.mutations.count()

@register(DrugRegimen)
class DrugRegimenAdmin(ModelAdmin):
    list_display = ('code', 'name', 'desc')

@register(DrugClass)
class DrugClassAdmin(ModelAdmin):
    list_display = ('code', 'name')

class MutationInline(StackedInline):
    model = Mutation
    extra = 2

@register(GeneLocus)
class GeneLocusAdmin(ModelAdmin):
    inlines = (MutationInline,)
    list_filter = ('genome', 'gene_type', 'strand')
    list_display = ('__str__', 'description', 'previous_id', 'gene_symbol', 'mutation_count', 'genome')
    search_fields = ('name', 'description', 'gene_symbol', 'previous_id')

    @staticmethod
    def mutation_count(obj):
        return obj.mutations.count()

site.register(Genome)

@register(Mutation)
class MutationAdmin(ModelAdmin):
    list_display = ('name', 'old_id', 'gene_locus', 'drugs_list')
    list_filter = ('drugs',)
    search_fields = ('name', 'old_id')

    @staticmethod
    def drugs_list(obj):
        return ", ".join(obj.drugs.values_list('code', flat=True))

site.register(TargetSet)
site.register(TargetRegion)

@register(StrainSource)
class StrainSourceAdmin(ModelAdmin):
    list_display = ('__str__', 'old_id', 'patient_id', 'country', 'date')
    list_filter = ('importer', 'source_lab', 'patient_sex', 'patient_hiv', 'resistance_group')
    search_fields = ('name', 'old_id', 'patient_id')

@register(StrainMutation)
class StrainMutationAdmin(ModelAdmin):
    search_fields = ('mutation__name', 'strain__name')
    list_display = ('mutation', 'strain')
    raw_id_fields = ('mutation', 'strain')

@register(StrainResistance)
class StrainResistanceAdmin(ModelAdmin):
    list_filter = ('resistance', 'drug')

site.register(BioProject)
site.register(Paper)
