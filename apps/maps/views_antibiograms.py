#
# Copyright (C) 2021  Dr. Maha Farhat
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
A second map with Marginal Resistance data.
"""

import re
import json
from collections import defaultdict
from django.views.generic import TemplateView
from django.db.models import Count

from apps.mutations.models import Drug

from .mixins import JsonView, DataSlicerMixin
from .utils import GraphData, adjust_coords
from .models import Country

from apps.maptables.models import CustomMap, MapRow

class AntibiogramMap(TemplateView):
    title = "Antibiograms Map"
    template_name = 'maps/map_antibiograms.html'

class MarginalPlaces(JsonView, DataSlicerMixin):
    """
    Provide a json output that slices mutation data by country
    """
    model = MapRow
    order = ['country__name', 'country__region']
    values = ['country__iso2', 'drug__code', 'drug__name', 'data']
    filters = {
        'parent_map__slug': ['antibiograms'],
        'drug[]': 'drug__code',
    }

    def get_context_data(self, **_):
        drugs = dict()
        rows = defaultdict(dict)
        for row in self.get_data():
            drugs[row['drug__code']] = row['drug__name']
            rows[row['country__iso2']][row['drug__code']] = json.loads(row['data'])

        # Limit to one drug only
        drug = list(drugs)[0]

        return {
            "type": "FeatureCollection",
            'filters': self.applied_filters(),
            'drugs': drugs,
            'features': [
                {
                    # Turning this to json and then back to python just to feed
                    # to JsonView, seems a little wasteful and redundent.
                    #"geometry": adjust_coords(json.loads(country.geom.geojson)),
                    "geometry": json.loads(country.geom.geojson),
                    "popupContent": country.name,
                    "type": "Feature",
                    "id": country.id,
                    "properties": {
                        "name": country.name,
                        "value": country.iso2,
                        "row": rows[country.iso2][drug],
                        "drug": {'name': drugs[drug], 'code': drug},
                    }

                } for country in Country.objects.filter(iso2__in=rows)
            ]
        }


class MarginalDrugs(JsonView, DataSlicerMixin):
    """Provide a json data slice into the drug resistance data"""
    model = MapRow
    order = ['drug__regimen', '-drug__priority',]
    values = ['drug__name', 'drug__code', 'data']
    filters = {
        'parent_map__slug': ['antibiograms'],
        'country[]': 'country__iso2__in',
    }

    def get_context_data(self, **_):
        """Return a dictionary of template variables"""

        totals = defaultdict(int)
        for row in self.get_data():
            dat = json.loads(row['data'])
            totals[row['drug__code']] += dat['gentb_snp10_n']

        
        drug_dict = [
            {
                'key': 'gentb_snp10_n',
                'values': [{
                    'y': i,
                    'x': t,
                    'col': t,
                    'value': i,
                    'total': i,
                } for t, i in totals.items()]
            }
        ]

        return {
            'data': drug_dict,
            'filters': self.applied_filters(),
        }

