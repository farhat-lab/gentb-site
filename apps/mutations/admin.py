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

