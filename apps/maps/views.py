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
from django.db.models import Count, Q, F, OuterRef, Subquery, Sum

from apps.mutations.models import (
    ImportSource, StrainSource, StrainMutation, GeneLocus, Genome, StrainResistance,
    Mutation, Paper, BioProject, Lineage, RESISTANCE, RESISTANCE_GROUP,
    StrainMutationIndex, StrainMutationCount
)

from .mixins import JsonView, DataSlicerMixin, DataTableMixin, PleaseWait
from .utils import GraphData, many_lookup, geo_adjust
from .models import Country, CountryHealth, CountryDetail

def get_gdp(self, **_):
    """ Getter for country detail gdp"""
    try:
        return self.detail.gdp
    except CountryDetail.DoesNotExist:
        return None

def get_health(self, **_):
    """ Getter for country health"""
    try:
        return self.health
    except CountryHealth.DoesNotExist:
        return CountryHealth()

def get_world_bank_gdp(self, **_):
    """ Getter for world bank gdp"""
    try:
        if not self.health.world_bank_gdp:
            return None
        return float(self.health.world_bank_gdp) / (10 ** 12)
    except CountryHealth.DoesNotExist:
        return None


class MapPage(TemplateView):
    """The html map page everything is provided by javascript"""
    title = "Mutations Map"
    template_name = 'maps/map_mutations.html'


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
        'lineage[]': 'lineage__name__in',
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
                    "geometry": geo_adjust(json.loads(country.geom.geojson)),
                    "popupContent": country.name,
                    "type": "Feature",
                    "id": country.id,
                    "properties": {
                        "name": country.name,
                        "value": country.iso2,
                        "values": ret[country.iso2],
                        "gdp": get_gdp(country),
                        "total_funding": get_health(country).total_funding,
                        "hiv_incidence2018": get_health(country).hiv_incidence2018,
                        "household": get_health(country).household,
                        "who_est_mdr": get_health(country).est_mdr,
                        "all_tb_incidence2018": get_health(country).all_tb_incidence2018,
                        "pop_dens": get_health(country).pop_dens,
                        "world_bank_gdp": get_world_bank_gdp(country),
                        "total_wealth": get_health(country).total_wealth
                    }

                } for country in Country.objects.filter(iso2__in=list(ret))
            ]
        }


class DrugList(JsonView, DataSlicerMixin):
    """Provide a json data slice into the drug resistance data"""
    model = StrainSource
    order = ['drugs__drug__regimen', '-drugs__drug__priority',]
    values = ['drugs__drug__name', 'drugs__drug__code', 'drugs__resistance']
    filters = {
        'lineage[]': 'lineage__name__in',
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
    """Get a list of loci that somewhat match the given locus string"""
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
    selected = ['mutation[]', 'name', str]

    cache_filters = {
        'source[]': 'importer_id',
        'paper[]': 'source_paper_id',
        'map[]': 'country__iso2',
        #'drug[]': 'drug__code',
        'lineage[]': 'lineage__name',
    }
    strain_filters = {
        'source[]': 'strain__importer__in',
        'paper[]': 'strain__source_paper__in',
        'map[]': 'strain__country__iso2__in',
        #'drug[]': 'strain__drug__code__in',
        'lineage[]': 'strain__lineage__name__in',
    }
    filters = {
        'genelocus[]': (int, 'gene_locus_id__in'),
    }

    def get_count_indexes(self):
        # Each combination of multiples (multiple selected items)
        # Must produce a cache for any of the caches to be valid
        indexes = []
        for query in self.filter_product(self.cache_filters):
            cache, created = StrainMutationIndex.objects.get_or_create(**query)
            #if not created:
            #    cache.requested += 1
            #    cache.save()
            indexes.append((cache, created))

        for cache, created in indexes:
            if cache.generating:
                raise PleaseWait("Generating mutations table, please check back soon...")
            elif created:
                raise PleaseWait("Requesting mutations table, please check back soon...")
            elif not cache.generated:
                raise PleaseWait("Mutations table requested, please check back soon...")

        return indexes

    def get_queryset(self):
        return super(Mutations, self).get_queryset().none()

    def hard_queryset(self, data):
        #qset = super(Mutations, self).get_queryset()

        query = self.apply_filters(self.strain_filters, prefix='strain_mutations__')
        qset = Mutation.objects.filter(query)
        qset = qset.annotate(strain_count=Count('strain_mutations__pk'),).filter(strain_count__gt=1)

        #query = self.apply_filters(self.strain_filters, Q(mutation=OuterRef('pk')))
        #strains = StrainMutation.objects.filter(query).order_by().values('mutation')
        #count_strains = strains.annotate(c=Count('*')).values('c')
        #qset = qset.annotate(strain_count=Subquery(count_strains)).filter(
        #    strain_count__gt=0,
        #)
        #qset = qset.annotate(strain_count=Count('strain_mutations__pk'),).filter(
        #    strain_count__gt=0,
        #)
        return qset

    def _process_datatable(self, qset, columns=(), order=(), search=None, start=0, length=-1, **_):
        """Special handling for the strain_count order_by case"""
        st_order = list(self.get_order(columns, order, strain_count='count', prefix='mutation__'))
        if 'count' in st_order or '-count' in st_order:
            # 1. If not searching, then we ignore text searches and use the index
            if search is None or not search.get('value', None):
                # XXX If the queryset is small, we could do the request directly.
                return self.process_counted_table(columns, int(start), int(length), st_order)

            # 2. If searching, then change the order, we don't want to search and order by count here.
            order = ()

        raise PleaseWait("Table is disabled for maintence")
        return super().process_datatable(qset,
            columns=columns, order=order, search=search, start=start, length=length)

    def process_counted_table(self, columns, start, length, order):
        """
        This special case for the strain order by allows us to use an index
        which gives us the pre-calculated counts.
        """
        # Resulting pile is a list of mutations
        db_columns = [self.column_to_django(col, prefix='mutation__', strain_count='count') for col in columns]
        counts = {}
        final_count = 0

        # We set a fixed number of pages to get, because they have to be processed every time.
        pages = length * 5

        for index, _ in self.get_count_indexes():
            # Get all the mutations, plus counts associated with this index, order first
            qset = index.mutations.order_by(*order)
            # Apply the same filters the process_datatable would have
            qset = qset.filter(self.apply_filters(self.filters, prefix='mutations__'))
            final_count += qset.count()

            # We get every item from zero to length, this is because of the combination
            for item in qset[0:pages].values('mutation_id', *db_columns):
                pk = item.pop('mutation_id')
                item = dict([(n.replace('mutation__', ''), v) for n ,v in item.items()])
                if pk in counts:
                    counts[pk]['strain_count'] += item['count']
                else:
                    item['strain_count'] = item.pop('count')
                    counts[pk] = item

        # Get selected items if found in this list of items
        selected_values = self.get_selected_values()
        selected = [counts.pop(pk) for pk in list(counts) if counts[pk]['name'] in selected_values]

        # Manually sort and slice the results
        ret = sorted(counts.values(),
            key=lambda item: item['strain_count'],
            reverse=('-count' in order)
        )
        end = start + length

        # This count is a fixed number of pages (5) or the length of the result
        # whichever is smaller.
        return selected, ret[start:end], min(len(ret), pages)

    def _prep_data(self, qset, columns, **extra):
        if qset is None:
            return []
        if isinstance(qset, list):
            return qset
        db_columns = [self.column_to_django(col) for col in columns]
        #if 'strain_count' in db_columns:
        #    db_columns.remove('strain_count')
        #    return [self.add_strain_count(item) for item in qset.values(*db_columns)]
        return list(qset.values(*db_columns))

    #def add_strain_count(self, row):
    #    """Manually count the strains in this row"""
    #    obj = Mutation.objects.get(name=row['name'])
    #    qset = obj.strain_mutations.filter(self.apply_filters(self.strain_filters))
    #    row['strain_count'] = qset.count()
    #    return row

    def post(self, request, *args, **kwargs):
        self.request.GET = self.request.POST
        return self.get(request, *args, **kwargs)


class MutationView(JsonView, DataSlicerMixin):
    """Provide a way to look at the resistance data via selected mutations"""
    model = StrainSource
    required = ['mutation[]']
    filters = {
        'source[]': 'importer__in',
        'paper[]': 'source_paper__in',
        'map[]': 'country__iso2__in',
    }

    @property
    def categories(self):
        """Return the categories available depending on the GET mode"""
        if 'drug[]' in self.request.GET:
            return OrderedDict(RESISTANCE)
        return OrderedDict(RESISTANCE_GROUP)

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""
        categories = dict([(str(key).upper(), value) for key, value in self.categories.items() if key])
        mutations = sorted(self.request.GET.getlist(self.required[0]))
        strains = self.get_data()

        (field, values) = many_lookup(StrainMutation, 'mutation__name', 'strain_id')(mutations)

        if 'drug[]' in self.request.GET:
            filters = self.applied_filters() + ['drug', 'mutation']
            columns = ['drugs__resistance', 'drugs__drug__code', 'mutations__mutation__name']

            drugs = sorted(self.request.GET.getlist('drug[]'))
            totals = []
            counts = []
            for drug in drugs:
                qs_totals = strains.filter(drugs__drug__code=drug)\
                                   .values_list(*columns[:-1])\
                                   .annotate(count=Count(columns[0]))

                qs_counts = strains.filter(drugs__drug__code=drug, **{field: values})\
                                   .values_list(*columns)\
                                   .annotate(count=Count(columns[0]))

                totals += [(str(row[0]).upper(), row[-1]) for row in qs_totals]
                # Multiple drugs get laid out (this is currently broken because
                # totals contains two of each category and so can't total up this
                # extra dimention correctly (fix in GraphData later)
                f_name = "{2} ({1})" if len(drugs) > 1 else "{2}"
                print(f"TOTALS for drug: {drug} is {totals}")

                counts += [{
                    'name': f_name.format(*row),
                    'count': row[-1],
                    'cat': str(row[0]).upper()
                } for row in qs_counts if row[-2] in mutations and row[1] == drug]

            mutations = sorted(list(set([row['name'] for row in counts])))
        else:
            filters = self.applied_filters() + ['mutation']
            columns = ['resistance_group', 'mutations__mutation__name']

            totals = strains.values_list(*columns[:-1])\
                            .annotate(count=Count(columns[0]))

            counts = strains.filter(**{field: values})\
                            .values_list(*columns)\
                            .annotate(count=Count(columns[0]))

            #orig = totals
            totals = [(str(row[0]).upper(), row[1]) for row in totals]
            counts = [{'name': row[1], 'count': row[2], 'cat': str(row[0]).upper()}
                      for row in counts if row[1] in mutations]

        return {
            'filters': filters,
            'data': GraphData(counts, 'name', 'count', 'cat', filter_label=self.filter_label)
                    .set_axis('z', categories, trim=True)
                    .set_axis('x', mutations, trim=False)
                    .set_axis('y', totals, trim=[None])
                    .to_graph()
            }

    @staticmethod
    def filter_label(axis, label):
        """Make sure filter labels are upper case"""
        if axis == 'z':
            return label.upper()
        return label
