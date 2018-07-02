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

import json
from collections import defaultdict, OrderedDict
from django.views.generic import TemplateView
from django.db.models import Count, Q, F, IntegerField
from django.db.models.functions import Cast

from apps.mutations.models import (
    ImportSource, StrainSource, Mutation, GeneLocus,
    RESISTANCE, RESISTANCE_GROUP,
)

from .mixins import JsonView, DataSlicerMixin
from .utils import GraphData
from .models import Country

LINEAGE_COLS = ['spoligotype_family', 'rflp_family', 'principle_group', 'wgs_group']
LINEAGE_NAMES = ['Spoligo', 'RFLP', 'PGG', 'WGS']

class MapPage(TemplateView):
    """The html map page everything is provided by javascript"""
    title = "Mutations Map"
    template_name = 'maps/map.html'


class Sources(JsonView, DataSlicerMixin):
    """
    Provide a json output that slices source import data
    """
    model = ImportSource
    order = ['pk']

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""
        _qs = self.get_data()
        return {
            'values': [
                [source.pk, source.name, str(source.uploader), source.strainsource_set.count()]
                for source in _qs
            ],
        }

class Places(JsonView, DataSlicerMixin):
    """
    Provide a json output that slices mutation data by country
    """
    model = StrainSource
    order = ['country__name', 'country__region']
    values = ['country__iso2', 'resistance_group']
    filters = dict(
        [
            ('drug', 'drugs__drug__code'),
            ('source', 'importer'),
        ] + zip(LINEAGE_NAMES, LINEAGE_COLS)
    )

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""
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
                    "properties": {
                        "name": country.name,
                        "value": country.iso2,
                        "values": ret[country.iso2],
                    },
                } for country in Country.objects.filter(iso2__in=list(ret))
            ],
        }


class Drugs(JsonView, DataSlicerMixin):
    """Provide a json data slice into the drug resistance data"""
    model = StrainSource
    order = ['drugs__drug__name', 'drugs__drug__kind']
    values = ['drugs__drug__name', 'drugs__drug__code', 'drugs__resistance']
    filters = dict(
        [
            ('map', 'country__iso2'),
            ('source', 'importer'),
        ] + zip(LINEAGE_NAMES, LINEAGE_COLS)
    )

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""
        return {
            'data': GraphData(
                self.get_data().annotate(count=Count('pk')),
                'drugs__drug__code', 'count', 'drugs__resistance',
            ).set_axis('z', RESISTANCE).to_graph(),
        }

class Lineages(JsonView, DataSlicerMixin):
    """Provide a json data slice into the Lineages data"""
    model = StrainSource
    order = ['spoligotype_family']
    values = LINEAGE_COLS
    filters = {
        'map': 'country__iso2',
        'drug': 'drugs__drug__code',
    }

    def get_queryset(self, without=None):
        _qs = super(Lineages, self).get_queryset()
        return _qs.filter(spoligotype_family__isnull=False)

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""
        return {
            'data': GraphData(
                self.get_data().annotate(count=Count('pk')),
                self.values, 'count', None, trim=True)
                    .set_axis('z', zip(self.values, LINEAGE_NAMES))
                    .set_axis('x', [(None, "Not Available")])
                    .to_graph()
        }

class Mutations(JsonView, DataSlicerMixin):
    """Provide a lookup into the mutations database for selecting anavailable mutation"""
    model = Mutation
    order = None
    values = ['pk']
    filters = {
        'snp': 'name__icontains',
        'ecoli': 'ecoli_aapos',
        'locus': 'gene_locus__name',
    }

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""
        if 'snp' not in self.request.GET and 'ecoli' not in self.request.GET:
            if 'locus' not in self.request.GET:
                return {'values': []}
            if 'range' in self.request.GET:
                return self.get_gene_range(**self.request.GET)
            return self.get_genes(**self.request.GET)

        # Otherwise mutation query
        qset = self.get_data()
        if qset.count() == 0:
            return {'msg': 'None found'}
        elif qset.count() > 200:
            return {'msg': 'Too many (%d)' % qset.count()}

        return {
            'msg': "Found %d mutations" % qset.count(),
            'values': list(self.get_my_list(qset)),
        }

    def get_gene_range(self, locus, synonymous=False, **_):
        """Returns a list of segments in a gene"""
        try:
            locus = GeneLocus.objects.get(name=locus[0])
        except GeneLocus.DoesNotExist:
            return {'title': "No Gene Locus found!"}

        mutations = Mutation.objects.filter(gene_locus=locus, nucleotide_position__isnull=False)
        if synonymous in (False, 0, 'false'):
            mutations = mutations.exclude(syn='S')
        qset = mutations.annotate(block=Cast(F('nucleotide_position')/50, IntegerField()))\
                        .values('block').annotate(count=Count('block'))
        offset = int(locus.start / 50.0)
        values = dict([(a - offset, b) for a, b in self.get_list(qset, 'block', 'count')])
        return {
            'start': locus.start,
            'end': locus.stop,
            'title': "{0.name} / {0.previous_id} ({1} mutations)".format(locus, mutations.count()),
            'values': values,
            'max': max(values.values()),
        }

    def get_genes(self, locus, **_):
        """Returns a list of genes"""
        locus = self.request.GET['locus']
        qset = GeneLocus.objects.filter(
            Q(name__istartswith=locus)
            | Q(previous_id__istartswith=locus)
            | Q(gene_symbol__istartswith=locus),
        )
        return {
            'msg': "Found %d genes" % qset.count(),
            'values': self.get_list(qset, 'name')
        }


    def get_my_list(self, _qs):
        """The core get list for thsi json data"""
        for (name, aar, eaa, aav) in self.get_list(
                _qs, 'name', 'aminoacid_reference', 'ecoli_aapos', 'aminoacid_varient'):
            if 'ecoli' in self.request.GET:
                yield {
                    'name': "%s+%s%s%s (E:%s)" % (name, aar, eaa, aav, eaa),
                    'value': name,
                }
            else:
                yield name



class MutationView(JsonView, DataSlicerMixin):
    """Provide a way to look at the resistance data via selected mutations"""
    model = StrainSource
    required = ['mutation[]',]
    filters = {
        'mutation[]': 'mutations__mutation__name__in',
        'drug': 'drugs__drug__code',
        'map': 'country__iso2',
        'src': 'importer',
    }
    @property
    def values(self):
        """Return drug or resistance values depending on the GET mode"""
        if 'drug' in self.request.GET:
            return ['mutations__mutation__name', 'drugs__resistance']
        return ['mutations__mutation__name', 'resistance_group']

    @property
    def categories(self):
        """Return the categories available depending on the GET mode"""
        if 'drug' in self.request.GET:
            return OrderedDict(RESISTANCE)
        return OrderedDict(RESISTANCE_GROUP)

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""
        mutations = self.request.GET.getlist(self.required[0])
        totals = self.get_data(without=self.values[0]).annotate(count=Count('pk'))
        totals = [(row[self.values[1]], row['count']) for row in totals]
        _qs = self.get_data().annotate(count=Count('pk'))
        return {
            'data': GraphData(_qs, self.values[0], 'count', self.values[-1])
                    .set_axis('z', self.categories, trim=True)
                    .set_axis('x', mutations)
                    .set_axis('y', totals, trim=[None])
                    .to_graph()
            }
