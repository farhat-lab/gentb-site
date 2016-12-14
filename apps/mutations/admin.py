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

class LocusInline(StackedInline):
    model = GeneLocus
    extra = 5

class DrugAdmin(ModelAdmin):
    inlines = (LocusInline,)

    def get_form(self, *args, **kw):
        return DrugForm

site.register(Drug, DrugAdmin)

class MutationInline(StackedInline):
    model = Mutation
    extra = 2

class GeneLocusAdmin(ModelAdmin):
    inlines = (MutationInline,)
    list_display = ('admin_name',)

    def admin_name(self, obj):
        return "%s for drug %s" % (obj.name, str(obj.drug))

site.register(GeneLocus, GeneLocusAdmin)

class MutationAdmin(ModelAdmin):
    list_display = ('name', 'gene_locus', 'drug')
    list_filter = ('gene_locus', 'gene_locus__drug')

    def drug(self, obj):
        return obj.gene_locus.drug

site.register(Mutation, MutationAdmin)

