#
# Copyright (C) 2017-2019 Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Test for maps
"""

import json
from collections import Counter

from django.test import TestCase
from extratest.base import ExtraTestCase

from ..mutations.models import ImportSource
from .utils import OrderlyDict, OrderedDict, GraphData

class UtilsTest(TestCase):
    """Test map utility functions"""
    def test_orderly(self):
        """Test OrderlyDict function"""
        item_a = OrderlyDict(OrderedDict([('a', 'a'), ('b', 'b'), ('c', 'c')]))
        self.assertEqual(item_a, OrderlyDict(['a', 'b', 'c']))

        item_b = OrderedDict([('a', 1), ('b', 2), ('c', 3)])
        self.assertEqual(item_b, OrderlyDict([('a', 1), ('b', 2), ('c', 3)]))

    def test_graph_basics(self):
        """Test GraphData function x and y only"""
        graph = list(GraphData([
            {'x': 'a', 'y': 2},
            {'x': 'a', 'y': 4},
            {'x': 'b', 'y': 10},
        ], 'x', 'y', None).to_graph())
        self.assertEqual(len(graph), 1) # Two z-axis categories
        values = graph[0]['values']
        self.assertEqual(len(values), 2) # Two y-axis values
        self.assertEqual(values[0], {'y': 6, 'x': 'a', 'total': -1, 'col': 'a', 'value': 6})
        self.assertEqual(values[1], {'y': 10, 'x': 'b', 'total': -1, 'col': 'b', 'value': 10})

    def test_graph_three_d(self):
        """Test GraphData function with z-axis"""
        graph = list(GraphData([
            {'x': 'a', 'y': 2, 'z': 'DOG'},
            {'x': 'a', 'y': 4, 'z': 'DOG'},
            {'x': 'b', 'y': 10, 'z': 'CAT'},
            {'x': 'a', 'y': 4, 'z': 'CAT'},
            {'x': 'b', 'y': 10, 'z': 'CAT'},
            {'x': 'b', 'y': 19, 'z': 'DUCK'},
            {'x': 'c', 'y': 42, 'z': 'DUCK'},
        ], 'x', 'y', 'z')\
            .set_axis('x', [('c', 'Carrot'), ('b', 'b-berry'), ('a', 'Apple')])\
            .set_axis('z', ['DOG', 'DUCK', 'CAT'])\
            .to_graph())

        self.assertEqual(len(graph), 3) # Three z-axis categories
        self.assertEqual(len(graph[0]['values']), 3) # Three x-axis values
        self.assertEqual(graph[0]['key'], 'DOG')
        self.assertEqual(len(graph[1]['values']), 3) # Three x-axis values
        self.assertEqual(graph[1]['key'], 'DUCK')
        self.assertEqual(len(graph[2]['values']), 3) # Three x-axis values (symetric)
        self.assertEqual(graph[2]['key'], 'CAT')

        values = graph[0]['values']
        self.assertEqual(values[2], {'y': 0, 'x': 'Carrot', 'total': -1, 'col': 'c', 'value': 0})
        self.assertEqual(values[0], {'y': 6, 'x': 'Apple', 'total': -1, 'col': 'a', 'value': 6})
        self.assertEqual(values[1], {'y': 0, 'x': 'b-berry', 'total': -1, 'col': 'b', 'value': 0})

        values = graph[2]['values']
        self.assertEqual(values[2], {'y': 0, 'x': 'Carrot', 'total': -1, 'col': 'c', 'value': 0})
        self.assertEqual(values[1], {'y': 4, 'x': 'Apple', 'total': -1, 'col': 'a', 'value': 4})
        self.assertEqual(values[0], {'y': 20, 'x': 'b-berry', 'total': -1, 'col': 'b', 'value': 20})

        values = graph[2]['values']
        self.assertEqual(values[2], {'y': 0, 'x': 'Carrot', 'total': -1, 'col': 'c', 'value': 0})
        self.assertEqual(values[1], {'y': 4, 'x': 'Apple', 'total': -1, 'col': 'a', 'value': 4})
        self.assertEqual(values[0], {'y': 20, 'x': 'b-berry', 'total': -1, 'col': 'b', 'value': 20})


class BaseCase(ExtraTestCase):
    """Basic functions for Json dats testing"""
    fixtures = ['test-genetics', 'test-maps', 'test-strains']
    maxDiff = 60000

    def assertJson(self, *args, **kwargs): # pylint: disable=invalid-name
        """
        Process a GET request back into context data from the JsonResponse
        """
        field = kwargs.pop('field', 'data')
        filters = kwargs.pop('filters', None)
        content = json.loads(self.assertGet(*args, **kwargs).content)
        if filters is not None:
            self.assertEqual(tuple(content['filters']), tuple(filters))
        return content.get(field)


class SourcesData(BaseCase):
    """Test sources data output (tab)."""
    def test_general_output(self):
        """Test output contains sources and papers only."""
        val = self.assertJson('maps:map.sources')['values']
        uni = Counter(["{kind}-{name}".format(**row) for row in val])
        self.assertEqual(len(uni), len(val), f"Sources aren't unique: {uni}")


class PlacesData(BaseCase):
    """Test places data output (tab)."""
    def assertPlaces(self, features, *tests): # pylint: disable=invalid-name
        """Check place data"""
        self.assertEqual(len(features), len(tests))
        for feature, (name, value, values) in zip(features, tests):
            self.assertEqual(feature['popupContent'], name)
            self.assertEqual(feature['properties']['name'], name)
            self.assertEqual(feature['properties']['value'], value)
            self.assertEqual(feature['properties']['values'], values)

    def test_all_output(self):
        """Test entire map output."""
        features = self.assertJson('maps:map.places', filters=(), field='features')

        self.assertTrue(features[0]['geometry']['coordinates'])
        self.assertEqual(features[0]['geometry']['type'], 'MultiPolygon')

        self.assertPlaces(
            features,
            ['France', 'FR', {'MDR': 3, 'Total': 3}],
            ['Germany', 'DE', {'MDR': 4, 'Total': 6, 'XDR': 1, 's': 1}],
            ['Russia', 'RU', {'MDR': 6, 'Total': 11, 'XDR': 4, 's': 1}],
        )

    def test_source_output(self):
        """Test source sliced map output."""
        source = ImportSource.objects.get(name='Import Z')
        maps = self.assertJson('maps:map.places', filters=('source',), field='features',
                               data={'source[]': [source.pk]})
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 1, 'Total': 1}],
            ['Germany', 'DE', {'MDR': 1, 'Total': 2, 'XDR': 1}],
            ['Russia', 'RU', {'MDR': 3, 'Total': 4, 's': 1}],
        )
        source = ImportSource.objects.get(name='Import Y')
        maps = self.assertJson('maps:map.places', filters=('source',), field='features',
                               data={'source[]': [source.pk]})
        self.assertPlaces(
            maps,
            ['Germany', 'DE', {'MDR': 3, 'Total': 4, 's': 1}],
            ['Russia', 'RU', {'MDR': 2, 'Total': 3, 'XDR': 1}],
        )
        source = ImportSource.objects.get(name='Import X')
        maps = self.assertJson('maps:map.places', filters=('source',), field='features',
                               data={'source[]': [source.pk]})
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 2, 'Total': 2}],
            ['Russia', 'RU', {'MDR': 1, 'Total': 4, 'XDR': 3}],
        )

    def test_paper_output(self):
        """Test paper sliced map output."""
        maps = self.assertJson('maps:map.places', filters=('paper',), data={'paper[]': [1]}, field='features')
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 2, 'Total': 2}],
            ['Germany', 'DE', {'MDR': 1, 'Total': 2, 'XDR': 1}],
            ['Russia', 'RU', {'MDR': 2, 'Total': 6, 'XDR': 4}],
        )
        maps = self.assertJson('maps:map.places', filters=('paper',), data={'paper[]': [2]}, field='features')
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 1, 'Total': 1}],
            ['Germany', 'DE', {'MDR': 3, 'Total': 4, 's': 1}],
            ['Russia', 'RU', {'MDR': 4, 'Total': 5, 's': 1}],
        )

    def test_drug_output_one(self):
        """Test drug sliced map output."""
        maps = self.assertJson('maps:map.places', filters=('drug',), data={'drug[]': ['H2O']}, field='features')
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 3, 'Total': 3}],
            ['Germany', 'DE', {'MDR': 4, 'Total': 5, 's': 1}],
            ['Russia', 'RU', {'MDR': 6, 'Total': 11, 'XDR': 4, 's': 1}],
        )

    def test_drug_output_many(self):
        """Test drug sliced map output."""
        maps = self.assertJson('maps:map.places', filters=('drug',), field='features',
                               data={'drug[]': ['MEM', 'WAVE', 'PIN', 'BUMP']})
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 3, 'Total': 3}],
            ['Germany', 'DE', {'MDR': 4, 'Total': 5, 's': 1}],
            ['Russia', 'RU', {'MDR': 6, 'Total': 11, 'XDR': 4, 's': 1}],
        )


class DrugListData(BaseCase):
    """Test drug data output."""
    def assertGraph(self, drugs, cols, rows): # pylint: disable=invalid-name
        """Test graph data"""
        for drug in drugs:
            self.assertEqual(rows[drug['key']],
                             [col['value'] for col in drug['values']])
            self.assertEqual(cols, [col['x'] for col in drug['values']])

    def test_everything(self):
        """Test unsliced drug output"""
        drugs = self.assertJson('maps:map.drugs', filters=())
        self.assertGraph(
            drugs,
            [None, 'BUMP', 'H2O', 'MEM', 'PIN', 'WAVE'], {
                'Unknown': [1, 0, 0, 0, 0, 0],
                'Sensitive to Drug': [0, 6, 9, 9, 8, 7],
                'Intermediate': [0, 0, 0, 0, 0, 0],
                'Resistant to Drug': [0, 13, 10, 9, 10, 11],
            })

    def test_source_output(self):
        """Test source sliced drug output"""
        source = ImportSource.objects.get(name='Import Z')
        drugs = self.assertJson('maps:map.drugs', filters=('source',),
                                data={'source[]': source.pk})
        self.assertGraph(
            drugs,
            [None, 'BUMP', 'H2O', 'MEM', 'PIN', 'WAVE'], {
                'Unknown': [1, 0, 0, 0, 0, 0],
                'Sensitive to Drug': [0, 3, 4, 5, 2, 3],
                'Intermediate': [0, 0, 0, 0, 0, 0],
                'Resistant to Drug': [0, 3, 2, 1, 4, 3],
            })

    def test_paper_output(self):
        """Test paper sliced drug output"""
        drugs = self.assertJson('maps:map.drugs', filters=('paper',), data={'paper[]': 1})
        self.assertGraph(
            drugs,
            [None, 'BUMP', 'H2O', 'MEM', 'PIN', 'WAVE'], {
                'Unknown': [1, 0, 0, 0, 0, 0],
                'Sensitive to Drug': [0, 2, 4, 3, 5, 2],
                'Intermediate': [0, 0, 0, 0, 0, 0],
                'Resistant to Drug': [0, 7, 5, 6, 4, 7],
            })

    def test_map_output(self):
        """Test map sliced drug output"""
        pass

class LineageData(BaseCase):
    """
    Test lineage data output.

     * Slice by source
     * Slice by paper
     * Slice by map country
     * Slice by drug
    """
    pass

class LocusRangeData(BaseCase):
    """
    Test locus range data output, should only show gene locuses with resistance data.

     * Slice by source
     * Slice by paper
     * Slice by map country
     * Slice by drug
    """
    pass

class LocusListData(BaseCase):
    """
    Test list of locus names.

     * Slice by Gene Type (Intergenic, Promoter, etc)
    """
    pass

class MutationsData(BaseCase):
    """
    Test mutations output, should only show mutations with resistances.

     * Slice by source
     * Slice by paper
     * Slice by map country
     * Slice by drug
     * Slice by gene locus
     * Slice by nucleotide position (range)
    """
    pass

class MutationResistanceData(BaseCase):
    """
    Test mutation resistance data.

     * Slice by source
     * Slice by paper
     * Slice by map country
     * Slice by drug
     * Slice by nucleotide position (range)
    """
    pass
