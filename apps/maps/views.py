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

import re
import json
from collections import defaultdict, OrderedDict
from django.views.generic import TemplateView, ListView
from django.views.decorators.csrf import csrf_exempt
from django.utils.decorators import method_decorator
from django.db.models import Count, Q

from apps.mutations.models import (
    ImportSource, StrainSource, Mutation, GeneLocus, Genome, StrainResistance,
    Paper, BioProject, Lineage, RESISTANCE, RESISTANCE_GROUP,
)

from .mixins import JsonView, DataSlicerMixin, DataTableMixin
from .utils import GraphData, many_lookup
from .models import Country, CountryHealth

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
        return {
            'values': list(self.get_sources()),
            'filters': self.applied_filters(),
        }

    def get_sources(self):
        """Return a list of data sources"""
        if self.request.GET.get('fields', '') == 'bio':
            for bioproject in BioProject.objects.filter(strains__isnull=False).distinct():
                yield dict(kind='bioproject', pk=bioproject.pk, name=bioproject.name,
                           count=bioproject.strains.count())
        else:
            for source in self.get_data():
                yield dict(kind='source', pk=source.pk, name=source.name,
                           uploader=str(source.uploader), count=source.strainsource_set.count())
            for paper in Paper.objects.filter(strains__isnull=False).distinct():
                yield dict(kind='paper', pk=paper.pk, name=paper.name,
                           url=paper.url, count=paper.strains.count())


class Places(JsonView, DataSlicerMixin):
    """
    Provide a json output that slices mutation data by country
    """
    model = StrainSource
    order = ['country__name', 'country__region']
    values = ['country__iso2', 'resistance_group']
    filters = {
        'source[]': 'importer__in',
        'paper[]': 'source_paper__in',
        'drug[]': many_lookup(StrainResistance, 'drug__code', 'strain_id'),
    }

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
            'filters': self.applied_filters(),
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


class DrugList(JsonView, DataSlicerMixin):
    """Provide a json data slice into the drug resistance data"""
    model = StrainSource
    order = ['drugs__drug__regimen', '-drugs__drug__priority',]
    values = ['drugs__drug__name', 'drugs__drug__code', 'drugs__resistance']
    filters = {
        'map[]': 'country__iso2__in',
        'source[]': 'importer__in',
        'paper[]': 'source_paper__in',
    }

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""

        drug_dict = list(GraphData(
            self.get_data().annotate(count=Count('pk')),
            'drugs__drug__code', 'count', 'drugs__resistance'
        ).set_axis('z', RESISTANCE).to_graph())

        # Sorting alphabetically by drug codename to prevent floating-bar errors in D3
        colors = ('#969696', '#3182bd', '#6baed6', '#fc8740')
        for x, section in enumerate(drug_dict):
            section['values'].sort(key=lambda el: (el['x'] or ''))
            # Set the colors manually
            section['color'] = colors[x]
            section['d'] = section['key'].lower().replace(' ', '_')

        if self.request.user.is_staff and \
                len(self.request.GET.getlist('map[]')) == 1:
            country = Country.objects.get(iso2=self.request.GET.get('map[]'))
            try:
                self.add_estimate_corrections(drug_dict, country.health.est_mdr / 100)
            except CountryHealth.DoesNotExist:
                section['COR'] = False

        return {
            'data': drug_dict,
            'filters': self.applied_filters(),
        }

    def add_estimate_corrections(self, graph, expected):
        """Correct for sensitivity error"""
        cors = []
        for x, drug in enumerate(graph[0]['values']):
            vals = dict([(section['d'], section['values'][x]['value']) for section in graph])
            # Calculate correction based on expected percentage.
            cor = self.estimate_correction(expected, **vals)
            cors.append({
                'x': drug['x'],
                'col': drug['col'],
                'value': cor,
                'y': cor,
                'total': -1,
            })
        graph.insert(1, {
            'values': cors,
            'color': "#9ecae1",
            'key': 'Oversampling',
            'expected': expected,
        })

    def estimate_correction(self, expected, sensitive_to_drug=0, intermediate=0, resistant_to_drug=0, **kw):
        """Estimate the corrects"""
        sensitive = sensitive_to_drug
        resistant = intermediate + resistant_to_drug
        total = float(sensitive + resistant)
        if total == 0 or not expected:
            return 0
        return max([int((resistant / expected) - total), 0])


class Lineages(JsonView, DataSlicerMixin):
    """
    Breakdown lineages with strain data added on.
    """
    model = StrainSource
    order = ['lineage__slug']
    values = ['lineage__name']
    filters = {
        'map[]': 'country__iso2__in',
        'source[]': 'importer__in',
        'paper[]': 'source_paper__in',
        'drug[]': many_lookup(StrainResistance, 'drug__code', 'strain_id'),
    }

    def get_sublineages(self, lin):
        """Returns list of all ancestors of `lin`, including `lin` (e.g. '3.6.4' -> ['3', '3.6', '3.6.4'])"""
        dot_indices = [idx for idx, ch in enumerate(lin) if ch == '.']
        return [lin[:idx] for idx in dot_indices] + [lin]

    def child_index(self, lin, children):
        """Returns index of child in `children` whose lineage is `lin`"""
        for idx, child in enumerate(children):
            if child['name'][1:] == lin:
                return idx


    def get_color(self, depth, idx):
        """Returns color of sunburst arc given depth and clockwise index.
           Arcs follow a blue->green gradient clockwise, with lower opacity at higher levels"""
        if depth == 0:
            return 'rgb(48,129,189)'
        blue_green = ['48,129,189', '49,147,185', '49,163,182', '49,178,179',
                      '49,175,158', '49,172,138', '49,169,119', '49,166,101', '48,163,84']
        return 'rgba({},{})'.format(blue_green[idx % 9], 1-depth/5.0)

    def lineage_tree(self):
        """
        Returns hierarchical dictionary of lineages for use in D3 charts.
        """

        # Contains frequency of each lineage (e.g. '4.3': 7)
        lin_counts = defaultdict(int)
        for strain in self.get_data():
            name = strain['lineage__name']
            if name and re.compile(r'[[A-Z0-9.]').match(name):
                lin_counts[name] += 1
        lin_counts = sorted(lin_counts.items())

        # Stores all lineages added to `lin_tree`
        processed = set()
        lin_tree = {
            'name': 'Total',
            'children': [],
            'filters': self.applied_filters(),
        }
        for lin, count in lin_counts:
            curr = lin_tree
            for depth, sl in enumerate(self.get_sublineages(lin)):
                if sl in processed:
                    next_idx = self.child_index(sl, curr['children'])
                    curr = curr['children'][next_idx]
                else:
                    curr['children'].append({'name': 'L'+sl, 'color': self.get_color(depth, len(curr['children'])), 'children': []})
                    curr = curr['children'][-1] # Step into last child
                    processed.add(sl)
            curr['size'] = count
        return lin_tree


    def get_context_data(self, **_):
        return self.lineage_tree()

class LocusList(DataTableMixin, ListView):
    """Get a list of locuses that somewhat match the given locus string"""
    model = GeneLocus
    search_fields = ['name', 'gene_symbol', 'description']
    filters = {}
    selected = ['genelocus[]', 'pk', int]

    def get_queryset(self):
        qset = super(LocusList, self).get_queryset()
        qset = qset.filter(start__isnull=False, stop__isnull=False)
        return qset.annotate(mcount=Count('mutations'))

@method_decorator(csrf_exempt, name='dispatch')
class Mutations(DataTableMixin, ListView):
    """Provide a lookup into the mutations database for selecting anavailable mutation"""
    model = Mutation
    search_fields = ['name', 'old_id', 'gene_locus__name']
    filters = {
        'source[]': 'strain_mutations__strain__importer__in',
        'paper[]': 'strain_mutations__strain__source_paper__in',
        'map[]': 'strain_mutations__strain__country__iso2__in',
        'drug[]': many_lookup(StrainResistance, 'drug__code', 'strain_id', 'strain_mutations__strain_id__in'),
        #'drug[]': 'strain_mutations__strain__drugs__drug__code__in',
        'genelocus[]': (int, 'gene_locus_id__in'),
    }
    selected = ['mutation[]', 'name', str]

    strain_filters = {
        'source[]': 'strain__importer__in',
        'paper[]': 'strain__source_paper__in',
        'drug[]': many_lookup(StrainResistance, 'drug__code', 'strain_id', 'strain_id__in'),
        'map[]': 'strain__country__iso2__in',
    }

    def prep_data(self, qset, columns, **kwargs):
        """Re-add counts for strains"""
        rows = super().prep_data(qset, columns, **kwargs)
        # Should only be 10 items as a maximum, meaning it should be 'ok' (but not great)
        for x, row in enumerate(rows):
            if x > 15:
                break
            obj = Mutation.objects.get(name=row['name'])
            # This is a manual smacking because the database can't count
            qset = obj.strain_mutations.filter(self.apply_filters(self.strain_filters))
            row['strain_count'] = qset.count()
        return rows

    def get_queryset(self):
        qset = super(Mutations, self).get_queryset()
        # We have to get the count so we can remove mutations with no strains too.
        qset = qset.annotate(strain_count=Count('strain_mutations__pk'),).filter(
            strain_count__gt=0,
        )
        return qset

    def post(self, request, *args, **kwargs):
        self.request.GET = self.request.POST
        return self.get(request, *args, **kwargs)

class MutationView(JsonView, DataSlicerMixin):
    """Provide a way to look at the resistance data via selected mutations"""
    model = StrainSource
    required = ['mutation[]',]
    filters = {
        'source[]': 'importer__in',
        'map[]': 'country__iso2__in',
        'drug[]': many_lookup(StrainResistance, 'drug__code', 'strain_id'),
        'mutation[]': 'mutations__mutation__name__in',
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
        totals = [(row[self.values[1]].upper(), row['count']) for row in totals]
        _qs = self.get_data().annotate(count=Count('pk'))
        return {
            'filters': self.applied_filters(),
            'data': GraphData(_qs, self.values[0], 'count', self.values[-1], filter_label=self.filter_label)
                    .set_axis('z', self.categories, trim=True)
                    .set_axis('x', mutations)
                    .set_axis('y', totals, trim=[None])
                    .to_graph()
            }

    def filter_label(self, axis, label):
        if axis == 'z':
            return label.upper()
        return label
