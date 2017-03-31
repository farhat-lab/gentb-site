from __future__ import print_function

import json
from collections import defaultdict, OrderedDict
from django.views.generic import TemplateView
from django.db.models import Count


from .json_view import JsonView
from .models import Country, Place
from apps.mutations.models import Drug, StrainSource, GeneLocus, Mutation, RESISTANCE

LINEAGE_COLS = ['spoligotype_family', 'rflp_family', 'principle_group', 'wgs_group']
LINEAGE_NAMES= ['Spoligo', 'RFLP', 'PGG', 'WGS']

class MapPage(TemplateView):
    template_name = 'maps/map.html'

class DataSlicerMixin(object):
    """
    Provide a way to slice up a given model based on inputs.
    """
    filters = {}
    order = []
    values = []

    def get_model(self):
        """Return the basic model to slice"""
        try:
            return self.model
        except:
            raise NotImplementedError("You must provide a model to slice.")

    def get_queryset(self):
        """Applies any filters from the request query to the given model"""
        def get_value(filtr, key):
            if filtr.endswith('__in'):
                return self.request.GET.getlist(key)
            return self.request.GET.get(key)

        qs = self.get_model().objects.all()
        filters = dict([(filtr, get_value(filtr, key))
            for (key, filtr) in self.filters.items()
                if key in self.request.GET])
        if self.order:
            qs = qs.order_by(*self.order)
        return qs.filter(**filters)

    def get_data(self):
        qs = self.get_queryset()
        if self.values:
            return qs.values(*self.values)
        return qs

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

class GraphData(list):
    """Format three columns into a format suitable for d3 graphs"""
    def __init__(self, qs, x, y, z, x_keys={}, z_keys={}):
        data = defaultdict(lambda: defaultdict(int))
        cols = OrderedDict()
        for dd in qs:
            # Collapse multiple fields into categories
            if isinstance(x, list):
                for tx in x:
                    data[tx][dd[tx]] += dd[y]
            # Or take categories from one field
            elif dd[y] > 0:
                cols[dd[x]] = 1
                data[dd.get(z, None)][dd[x]] += dd[y]
        
        # Make the data structure square and convert from defaultdicts to OrderedDicts
        for key in data:
            ret2 = []
            for col in (cols or data[key]):
                ret2.append({
                  "x": x_keys.get(col, col),
                  "y": data[key][col],
                })
            self.append({
              "key": z_keys.get(key, key),
              "values": ret2,
            })


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
            z_keys=dict(RESISTANCE),
          )
        }

class Lineages(JsonView, DataSlicerMixin):
    model = StrainSource
    order = ['spoligotype_family']
    values = LINEAGE_COLS
    filters = {
      'map': 'country__iso2',
      'drug': 'drugs__drug__code',
    }

    def get_queryset(self):
        qs = super(Lineages, self).get_queryset()
        return qs.filter(spoligotype_family__isnull=False)

    def get_context_data(self, **kw):
        return {
          'data': GraphData(
            self.get_data().annotate(count=Count('pk')),
            self.values, 'count', True,
            z_keys=dict(zip(self.values, LINEAGE_NAMES)),
            x_keys={None: "Not Available"},
          )
        }

class Mutations(JsonView, DataSlicerMixin):
    model = Mutation
    order = None
    values = ['pk']
    filters = {
      'drug': 'drugs__code',
    }

    def get_context_data(self, **kw):
        qs = self.get_data()

        ret = { 
          'levels': ['Gene Locus', 'Mutation'],
          'children': [], 
        }

        mutations = qs[:1000].values_list('name', 'gene_locus__name')

        out = defaultdict(list)
        for mutation, locus in mutations:
            out[locus].append(mutation)

        for locus in sorted(out):
            ret['children'].append({
              'name': locus,
              'children': [{'name': m, 'value': m} for m in out[locus]],
            })

        return ret 


        return {
          'data': GraphData(
            self.get_data()\
                .annotate(count=Count('strain_mutations__strain__pk'))
                [:20],
            'name', 'count', None,
            z_keys={None: 'All Mutations'},
            x_keys={None: "Not Available"},
          )
        }

class MutationView(JsonView, DataSlicerMixin):
    model = StrainSource
    values = ['mutations__mutation__name', 'resistance_group']
    filters = {
        'mutations[]': 'mutations__mutation__name__in',
            #  'map': 'country__iso2',
            #  'drug': 'drugs__drug__code',
    }

    def get_context_data(self, **kw):
        qs = self.get_data()
        if 'mutations[]' not in self.request.GET:
            qs = qs.filter(name='NOOP')
        #drug = self.request.GET.get('drug', None)
        #Mutation.objects.filter(name__in=self.snps).values('code')
        return {
           'data': GraphData(
               qs.annotate(count=Count('pk')),
                   'mutations__mutation__name', 'count', 'resistance_group',
                 )
               }

