#
# Copyright (C) 2017 Maha Farhat
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

from django.test import TestCase
from extratest.base import ExtraTestCase

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



class SourcesData(TestCase):
    """
    Test sources data output.

     * No specific slicing, full list only.

    """
    def test_general_output(self):
        """
        Test output contains sources and papers only.
        """
        pass

    def test_bioproject_output(self):
        """
        Test listing BioProject sources
        """
        pass

class PlacesData(TestCase):
    """
    Test places data output.

     * Slice by source
     * Slice by paper
     * Slice by drugs

    """
    def test_all_output(self):
        """
        Test entire map output.
        """
        pass

    def test_source_output(self):
        """
        Test source sliced map output.
        """
        pass

    def test_paper_output(self):
        """
        Test paper sliced map output.
        """
        pass

class DrugListData(TestCase):
    """
    Test drug data output.

     * Slice by source
     * Slice by paper
     * Slice by map country

    """
    pass

class LineageData(TestCase):
    """
    Test lineage data output.

     * Slice by source
     * Slice by paper
     * Slice by map country
     * Slice by drug
    """
    pass

class LocusRangeData(TestCase):
    """
    Test locus range data output, should only show gene locuses with resistance data.

     * Slice by source
     * Slice by paper
     * Slice by map country
     * Slice by drug
    """
    pass

class LocusListData(TestCase):
    """
    Test list of locus names.

     * Slice by Gene Type (Intergenic, Promoter, etc)
    """
    pass

class MutationsData(TestCase):
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

class MutationResistanceData(TestCase):
    """
    Test mutation resistance data.

     * Slice by source
     * Slice by paper
     * Slice by map country
     * Slice by drug
     * Slice by nucleotide position (range)
    """
    pass
