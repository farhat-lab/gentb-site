from __future__ import print_function

import json
from collections import defaultdict, OrderedDict
from django.views.generic import TemplateView
from django.db.models import Count


from .json_view import JsonView
from .models import Country, Place
from apps.mutations.models import Drug, StrainSource, RESISTANCE


class MapPage(TemplateView):
    template_name = 'maps/basic_map.html'

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
        qs = self.get_model().objects.all()
        filters = dict([(filtr, self.request.GET[key])
            for (key, filtr) in self.filters.items()
                if key in self.request.GET])
        return qs.filter(**filters)

    def get_data(self):
        qs = self.get_queryset()
        if self.order:
            qs = qs.order_by(*self.order)
        if self.values:
            return qs.values(*self.values)
        return qs

class Places(JsonView, DataSlicerMixin):
    model = Country
    order = ['name', 'region']
    values = ['iso2', 'sources__resistance_group']

    def get_context_data(self, **kw):
        ret = defaultdict(lambda: defaultdict(int))
        for row in self.get_data().annotate(count=Count('sources__pk')):
            group = row['sources__resistance_group']
            if group == 'S':
                group = 'Sensitive'
            if group is not None:
                ret[row['iso2']][group] = row['count']
                ret[row['iso2']]['Total'] += row['count']

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
              "properties": {"name": country.name, "values": ret[country.iso2]},
            } for country in Country.objects.filter(iso2__in=list(ret))
           ],
        }            


class Drugs(JsonView, DataSlicerMixin):
    model = Drug
    order = ['name', 'kind']
    values = ['name', 'code', 'strains__resistance']
    dr_key = dict(RESISTANCE)

    def get_context_data(self, **kw):
        data = defaultdict(lambda: defaultdict(int))
        cols = set()
        for dd in self.get_data().annotate(count=Count('strains__pk')):
            if dd['count'] > 0:
                cols.add(dd['code'])
                data[dd['strains__resistance']][dd['code']] += dd['count']
        
        # Make the data structure square and convert from defaultdicts to OrderedDicts
        ret = []
        for key in data:
            ret2 = []
            for col in cols:
                ret2.append({
                  "x": col,
                  "y": data[key][col],
                })
            ret.append({
              "key": self.dr_key[key],
              "values": ret2,
            })

        return {
          'drugs': ret,
        }

class Lineages(JsonView, DataSlicerMixin):
    model = StrainSource
    order = ['name', 'kind']
    values = ['spoligotype_family']

class Mutations(JsonView, DataSlicerMixin):
    pass

