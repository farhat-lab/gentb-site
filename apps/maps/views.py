#
# Copyright (C) 2017  Dr. Maha Farhat
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
Views for the mapping application
"""

from __future__ import print_function

import sys
import json
from collections import defaultdict, OrderedDict
from django.views.generic import TemplateView
from django.db.models import Count

from .mixins import JsonView, DataSlicerMixin
from .utils import GraphData
from .models import Country, Place
from apps.mutations.models import Drug, StrainSource, GeneLocus, Mutation, RESISTANCE, RESISTANCE_GROUP

LINEAGE_COLS = ['spoligotype_family', 'rflp_family', 'principle_group', 'wgs_group']
LINEAGE_NAMES= ['Spoligo', 'RFLP', 'PGG', 'WGS']

class MapPage(TemplateView):
    template_name = 'maps/map.html'

class Places(JsonView, DataSlicerMixin):
    model = StrainSource
    order = ['country__name', 'country__region']
    values = ['country__iso2', 'resistance_group']
    filters = dict(
        [('drug', 'drugs__drug__code')] +
        zip(LINEAGE_NAMES, LINEAGE_COLS)
    )

    def get_context_data(self, **kw):
        ret = defaultdict(lambda: defaultdict(int))
        for row in self.get_data().annotate(count=Count('pk')):
            group = row['resistance_group']
            if group == 'S':
                group = 'Sensitive'
            if group is not None:
                ret[row['country__iso2']][group] = row['count']
                ret[row['country__iso2']]['Total'] += row['count']

        return {
          "type": "FeatureCollection",
          'features': [
            {
              # Turning this to json and then back to python just to feed
              # to JsonView, seems a little wasteful and redundent.
              "geometry": json.loads(country.geom.geojson),
              "popupContent": country.name,
              "type": "Feature",
              "id": country.id,
              "properties": {"name": country.name, "value": country.iso2, "values": ret[country.iso2]},
            } for country in Country.objects.filter(iso2__in=list(ret))
           ],
        }


class Drugs(JsonView, DataSlicerMixin):
    model = StrainSource
    order = ['drugs__drug__name', 'drugs__drug__kind']
    values = ['drugs__drug__name', 'drugs__drug__code', 'drugs__resistance']
    filters = dict(
        [('map', 'country__iso2')] +
        zip(LINEAGE_NAMES, LINEAGE_COLS)
    )

    def get_context_data(self, **kw):
        return {
          'data': GraphData(
            self.get_data().annotate(count=Count('pk')),
            'drugs__drug__code', 'count', 'drugs__resistance',
          ).set_axis('z', RESISTANCE).to_graph(),
        }

class Lineages(JsonView, DataSlicerMixin):
    model = StrainSource
    order = ['spoligotype_family']
    values = LINEAGE_COLS
    filters = {
      'map': 'country__iso2',
      'drug': 'drugs__drug__code',
    }

    def get_queryset(self, without=None):
        qs = super(Lineages, self).get_queryset()
        return qs.filter(spoligotype_family__isnull=False)

    def get_context_data(self, **kw):
        return {
          'data': GraphData(
            self.get_data().annotate(count=Count('pk')),
            self.values, 'count', None, trim=True)
              .set_axis('z', zip(self.values, LINEAGE_NAMES))
              .set_axis('x', [(None, "Not Available")])
              .to_graph()
        }

class Mutations(JsonView, DataSlicerMixin):
    model = Mutation
    order = None
    values = ['pk']
    filters = {
      'snp': 'name__icontains',
      'locus': 'gene_locus__name',
    }

    def get_context_data(self, **kw):
        qs = self.get_data()
        if 'snp' not in self.request.GET:
            return {
                'values': self.get_list(qs, 'gene_locus__name'),
            }
        elif qs.count() == 0:
            return {'msg': 'None found'}
        elif qs.count() > 200:
            return {'msg': 'Too many (%d)' % qs.count()}

        return {
            'msg': "Found %d mutations" % qs.count(),
            'values': self.get_list(qs, 'name'),
        }


class MutationView(JsonView, DataSlicerMixin):
    model = StrainSource
    required = ['mutation[]',]
    filters = {
        'mutation[]': 'mutations__mutation__name__in',
        'drug': 'drugs__drug__code',
        'map': 'country__iso2',
    }
    @property
    def values(self):
        if 'drug' in self.request.GET:
            return ['mutations__mutation__name', 'drugs__resistance']
        return ['mutations__mutation__name', 'resistance_group']

    @property
    def categories(self):
        if 'drug' in self.request.GET:
            return OrderedDict(RESISTANCE)
        return OrderedDict(RESISTANCE_GROUP)

    def get_context_data(self, **kw):
        mutations = self.request.GET.getlist(self.required[0])
        totals = self.get_data(without=self.values[0]).annotate(count=Count('pk'))
        totals = [(row[self.values[1]], row['count']) for row in totals]
        qs = self.get_data().annotate(count=Count('pk'))
        return {
          'data': GraphData(qs, self.values[0], 'count', self.values[-1])
            .set_axis('z', self.categories, trim=True)
            .set_axis('x', mutations)
            .set_axis('y', totals)
            .to_graph()
        }

