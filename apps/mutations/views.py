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
# pylint: disable=too-few-public-methods
#
"""
Provides the drug, locus, mutation data for the progressive dropdown inputs.
"""

from collections import defaultdict

from apps.maps.mixins import JsonView

from django.views.generic import TemplateView, FormView, DetailView, ListView

from .models import ImportSource, Mutation
from .forms import DataUploaderForm
from .utils import info_mutation_format, match_snp_name


class MutationView(TemplateView):
    """View mutations"""
    template_name = "mutations/mutation.html"

    def get_context_data(self, **kw):
        """Gather together info of mutations"""
        dat = super(MutationView, self).get_context_data(**kw)
        if 'mutation' in self.request.GET:
            dat['highlight'], dat['regular_exp'], dat['info'] =\
                info_mutation_format(self.request.GET['mutation'])
        return dat


class UploadData(FormView):
    """Upload a new data importer."""
    title = "Upload Data to GenTB Mutations Tracker"
    template_name = 'mutations/upload_data.html'
    form_class = DataUploaderForm

    def get_success_url(self):
        """Return the Status page for the uploaded data"""
        return self.object.get_absolute_url()

    def form_valid(self, form):
        """The form isvalid so we'll save our form"""
        self.object = form.save(self.request.user) # pylint: disable=attribute-defined-outside-init
        return super(UploadData, self).form_valid(form)


class UploadView(DetailView):
    """A view of uploads"""
    template_name = 'mutations/upload_status.html'
    model = ImportSource


class UploadList(ListView):
    """A way to list a uploaded imports"""
    model = ImportSource
    title = "My Imports"
    parent = ("/maps/", "Maps")

    def get_queryset(self):
        """Return all import sources owned by the logged in user"""
        qset = super(UploadList, self).get_queryset()
        if not self.request.user.is_staff:
            qset = qset.filter(uploader=self.request.user)
        return qset

class DropDownData(JsonView):
    """Drop and drag data response"""
    @staticmethod
    def get_context_data():
        """Gather together all the data for drugs"""
        ret = {
            'levels': ['Gene Locus', 'Mutation'],
            'children': [],
        }
        names = Mutation.objects.all().variant_names()
        matrix = defaultdict(list)
        for name in names:
            lookup = match_snp_name(name)
            if 'rgene' in lookup and lookup['rgene']:
                lookup['gene'] = lookup['rgene']
            if 'noncode' in lookup and lookup['noncode']:
                if not lookup['noncode'].isdigit():
                    lookup['gene'] = lookup['noncode'] + ' ' + lookup['gene']
            if 'gene' in lookup and lookup['gene']:
                matrix[lookup['gene']].append(name)
            else:
                raise ValueError("Gene name missing from: {}".format(name))

        for gene in sorted(matrix):
            ret['children'].append({
                'name': "{}".format(str(gene)),
                'children': [{'name': name, 'value': name} for name in sorted(matrix[gene])],
            })
        return ret
