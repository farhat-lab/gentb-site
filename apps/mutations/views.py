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
Provides the drug, locus, mutation data for the progressive dropdown inputs.
"""

import logging
LOGGER = logging.getLogger(__name__)

from django.views.generic import *
from django.core.urlresolvers import reverse

from apps.maps.mixins import JsonView

from .models import *
from .forms import DataUploaderForm
from .utils import info_mutation_format

class MutationView(TemplateView):
    template_name = "mutations/mutation.html"

    def get_context_data(self, **kw):
        d = super(MutationView, self).get_context_data(**kw)
        if 'mutation' in self.request.GET:
            (hl, rs, info) = info_mutation_format(self.request.GET['mutation'])
            d['highlight'] = hl
            d['regular_exp'] = rs
            d['info'] = info
        return d


class UploadData(FormView):
    title = "Upload Data to GenTB Mutations Tracker"
    #parent = ("genes:upload.list", "My Imports")
    template_name = 'mutations/upload_data.html'
    form_class = DataUploaderForm

    def get_success_url(self):
        return self.object.get_absolute_url()

    def form_valid(self, form):
        self.object = form.save(self.request.user)
        return super(UploadData, self).form_valid(form)


class UploadView(DetailView):
    #parent = ("genes:upload.list", "My Imports")
    template_name = 'mutations/upload_status.html'
    model = ImportSource


class UploadList(ListView):
    model = ImportSource
    title = "My Imports"
    parent = ("/maps/", "Maps")

    def get_queryset(self):
        qs = super(UploadList, self).get_queryset()
        return qs.filter(uploader=self.request.user)


class DropDownData(JsonView):
    def get_context_data(self, *kw):
        ret = {
          'levels': ['Drug', 'Gene Locus', 'Mutation'],
          'children': [],
        }
        for drug in Drug.objects.all():
            ret['children'].append({
                'name': str(drug),
                'children': [],
            })
            qs = drug.mutations.all()
            if not self.request.GET.get('all', False):
                qs = qs.filter(predictor=True)
            loci = qs.values_list('gene_locus__name', flat=True).distinct()
            for locus in GeneLocus.objects.filter(name__in=loci):
                ret['children'][-1]['children'].append({
                  'name': str(locus),
                  'children': [],
                })
                for mutation in qs.filter(gene_locus__name=locus):
                    ret['children'][-1]['children'][-1]['children'].append({
                      'name': str(mutation),
                      'value': mutation.name,
                    })
        return ret

