from __future__ import print_function

from django.views.generic import TemplateView

import json

from .json_view import JsonView
from .models import Country, Place

class MapPage(TemplateView):
    template_name = 'maps/basic_map.html'

class Countries(JsonView):
    def get_context_data(self, **kw):
        return {
          "type": "FeatureCollection",
          'features': [
            {
              # Turning this to json and then back to python just to feed
              # to JsonView, seems a little wasteful and redundent.
              "geometry": json.loads(country.geom.geojson),
              "popupContent": country.name,
              "type": "Feature",
              "id": country.iso3,
              "properties": {"name": country.name},
            } for country in Country.objects.exclude(sources__isnull=True)
           ],
        }            

