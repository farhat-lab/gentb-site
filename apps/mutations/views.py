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

from apps.maps.json_view import JsonView

from .models import *

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

