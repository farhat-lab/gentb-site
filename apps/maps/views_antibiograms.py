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
from django.views.generic import ListView, TemplateView
from django.db.models import Count

from apps.mutations.models import Drug

from .mixins import JsonView, DataSlicerMixin
from .utils import GraphData, geo_adjust
from .models import Country

from apps.maptables.models import MapDisplay, MapDataRow
from apps.mutations.models import Drug

class AntibiogramMap(ListView):
    model = MapDisplay
    title = "Antibiograms Map"
    template_name = 'maps/map_antibiograms.html'

class MarginalPlaces(JsonView, DataSlicerMixin):
    """
    Provide a json output that slices mutation data by country
    """
    model = MapDataRow
    order = ['country__name', 'country__region']
    values = ['country__iso2', 'drug__code', 'drug__name', 'data']
    filters = {
        'drug[]': 'drug__code',
    }

    def get_context_data(self, **_):
        drugs = dict()
        parent_map = MapDisplay.objects.get(slug=self.kwargs['slug'])
        filters = list(self.prepare_filters(parent_map.data_filters.all()))
        rows = defaultdict(dict)
        values = []
        for row in self.get_data().filter(source_id=parent_map.data.pk):
            drugs[row['drug__code']] = row['drug__name']
            data = json.loads(row['data'])
            if self.filter_row(filters, data):
                rows[row['country__iso2']][row['drug__code']] = data
                if parent_map.fill_column in data:
                    try:
                        values.append(float(data[parent_map.fill_column]))
                    except Exception:
                        pass

        drug = None
        default_drug = 'INH'
        if parent_map.default_drug:
            default_drug = parent_map.default_drug.code
        drug_sort = [default_drug]

        # Limit to one drug only
        if len(drugs) > 1:
            for dname in drug_sort:
                if dname in drugs:
                    drug = dname
            if drug is None:
                drug = list(drugs)[0]
        elif len(drugs) == 1:
            drug = list(drugs)[0]

        all_drugs = dict(parent_map.data.rows.values_list('drug__code', 'drug__name'))
        fill_max = parent_map.get_max(values)

        return {
            "type": "FeatureCollection",
            "fill": {
                'max': fill_max,
                'column': parent_map.fill_column,
                'ranges': parent_map.get_ranges(fill_max),
                'colors': parent_map.get_colors(),
            },
            "details": [
                {'label': item.label, 'column': item.column, 'type': item.kind}
                    for item in parent_map.details.all()
            ],
            'c_filters': self.drug_filter(drug, all_drugs) + list(self.applied_filters(filters)),
            'features': [
                {
                    "srid": country.geom.srid,
                    "geometry": geo_adjust(json.loads(country.geom.geojson)),
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

    def drug_filter(self, drug, drugs):
        """Drug filter always available"""
        return [
            {'label': 'Drug', 'default': drug, 'key': 'drug', 'values': 
                [{'value': code, 'label': f"{label} ({code})", 'selected': (code == drug)}
                    for code, label in drugs.items()]
            },
        ]

    def prepare_filters(self, filters):
        """Prepare filter functions and values"""
        for _fl in filters:
            value = self.request.GET.getlist(_fl.key + '[]', self.request.GET.get(_fl.key, None))
            yield (_fl, getattr(self, 'filter_' + _fl.kind, self.filter_none), value)

    def filter_row(self, filters, row):
        for _fl, func, value in filters:
            if value is None or len(value) == 0:
                return True
            if not func(_fl, value[0], row):
                return False
        return True

    def filter_none(self, _fl, value, row):
        return True

    def filter_limit(self, _fl, value, row):
        """Filter based on the row"""
        if str(value).isdigit():
            value = int(value)
        if value is None:
            return True
        if _fl.column not in row:
            return False
        if row[_fl.column] < value:
            return False
        return True

    def applied_filters(self, filters):
        """Generate filters which the javascript can create into dropdowns"""
        for _fl, _, selected in filters:
            opts = _fl.get_options()
            if not opts or isinstance(opts, dict) and opts['error']:
                print("Ignoring missing or error drop down")
                continue

            if selected is None or len(selected) == 0:
                selected = opts[0].get('value', None)

            if not isinstance(selected, (tuple, list)):
                selected = (selected,)

            for op in opts:
                op['selected'] = op.get('value', None) in selected

            yield {
                'key': _fl.key,
                'label': _fl.label,
                'values': opts,
            }
