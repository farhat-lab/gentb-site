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

from ..mutations.models import ImportSource, GeneLocus
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

    def assertGraph(self, got_rows, cols, rows): # pylint: disable=invalid-name
        """Test graph data"""
        for x, got_row in enumerate(got_rows):
            with self.subTest(row=x, key=got_row['key']):
                self.assertEqual(set(cols), set([col['x'] for col in got_row['values']]))
                by_col = dict([(row['col'], row) for row in got_row['values']])
                self.assertEqual(rows[got_row['key']],
                                 [by_col[col]['value'] for col in cols])

    def assertGraphTotals(self, got_rows, cols, rows): # pylint: disable=invalid-name
        """Test graph data with totals"""
        for x, got_row in enumerate(got_rows):
            with self.subTest(row=x, key=got_row['key']):
                self.assertEqual(set(cols), set([col['x'] for col in got_row['values']]))
                by_col = dict([(row['col'], row) for row in got_row['values']])
                self.assertEqual(rows[got_row['key']][0],
                                 [by_col[col]['value'] for col in cols])
                self.assertEqual(rows[got_row['key']][1],
                                 [by_col[col]['total'] for col in cols])

    def assertDataTable(self, url, names=('pk', 'str'), start=0, length=5, **kwargs): # pylint: disable=invalid-name
        """Test data table output"""
        data = kwargs.pop('data', {})
        data.update({
            'draw': 1,
            'search[value]': '',
            'order[0][column]': kwargs.pop('order', 0),
            'order[0][dir]': 'asc',
            'start': start,
            'length': length,

        })
        for x, name in enumerate(names):
            data[f'columns[{x}][data]'] = name

        return self.assertJson(url, data=data, **kwargs)

class SourcesData(BaseCase):
    """Test sources data output (tab)."""
    def test_general_output(self):
        """Test output contains sources and papers only."""
        val = self.assertJson('maps:data.sources', field='values')
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
        features = self.assertJson('maps:data.places', filters=(), field='features')

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
        maps = self.assertJson('maps:data.places', filters=('source',), field='features',
                               data={'source[]': [source.pk]})
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 1, 'Total': 1}],
            ['Germany', 'DE', {'MDR': 1, 'Total': 2, 'XDR': 1}],
            ['Russia', 'RU', {'MDR': 3, 'Total': 4, 's': 1}],
        )
        source = ImportSource.objects.get(name='Import Y')
        maps = self.assertJson('maps:data.places', filters=('source',), field='features',
                               data={'source[]': [source.pk]})
        self.assertPlaces(
            maps,
            ['Germany', 'DE', {'MDR': 3, 'Total': 4, 's': 1}],
            ['Russia', 'RU', {'MDR': 2, 'Total': 3, 'XDR': 1}],
        )
        source = ImportSource.objects.get(name='Import X')
        maps = self.assertJson('maps:data.places', filters=('source',), field='features',
                               data={'source[]': [source.pk]})
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 2, 'Total': 2}],
            ['Russia', 'RU', {'MDR': 1, 'Total': 4, 'XDR': 3}],
        )

    def test_paper_output(self):
        """Test paper sliced map output."""
        maps = self.assertJson('maps:data.places', filters=('paper',), data={'paper[]': [1]}, field='features')
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 2, 'Total': 2}],
            ['Germany', 'DE', {'MDR': 1, 'Total': 2, 'XDR': 1}],
            ['Russia', 'RU', {'MDR': 2, 'Total': 6, 'XDR': 4}],
        )
        maps = self.assertJson('maps:data.places', filters=('paper',), data={'paper[]': [2]}, field='features')
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 1, 'Total': 1}],
            ['Germany', 'DE', {'MDR': 3, 'Total': 4, 's': 1}],
            ['Russia', 'RU', {'MDR': 4, 'Total': 5, 's': 1}],
        )

    def test_drug_output_one(self):
        """Test drug sliced map output."""
        maps = self.assertJson('maps:data.places', filters=('drug',), data={'drug[]': ['H2O']}, field='features')
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 3, 'Total': 3}],
            ['Germany', 'DE', {'MDR': 4, 'Total': 5, 's': 1}],
            ['Russia', 'RU', {'MDR': 6, 'Total': 11, 'XDR': 4, 's': 1}],
        )

    def test_drug_output_many(self):
        """Test drug sliced map output."""
        maps = self.assertJson('maps:data.places', filters=('drug',), field='features',
                               data={'drug[]': ['MEM', 'WAVE', 'PIN', 'BUMP']})
        self.assertPlaces(
            maps,
            ['France', 'FR', {'MDR': 3, 'Total': 3}],
            ['Germany', 'DE', {'MDR': 4, 'Total': 5, 's': 1}],
            ['Russia', 'RU', {'MDR': 6, 'Total': 11, 'XDR': 4, 's': 1}],
        )


class DrugListData(BaseCase):
    """Test drug data output."""
    def test_everything(self):
        """Test unsliced drug output"""
        drugs = self.assertJson('maps:data.drugs', filters=())
        self.assertGraph(
            drugs,
            [None, 'BUMP', 'H2O', 'MEM', 'PIN', 'WAVE'], {
                'Unknown': [1, 0, 0, 0, 0, 0],
                'Sensitive to Drug': [0, 6, 9, 9, 8, 7],
                'Intermediate': [0, 0, 0, 0, 0, 0],
                'Resistant to Drug': [0, 13, 10, 9, 8, 11],
            })

    def test_source_output(self):
        """Test source sliced drug output"""
        source = ImportSource.objects.get(name='Import Z')
        drugs = self.assertJson('maps:data.drugs', filters=('source',),
                                data={'source[]': source.pk})
        self.assertGraph(
            drugs,
            [None, 'BUMP', 'H2O', 'MEM', 'PIN', 'WAVE'], {
                'Unknown': [1, 0, 0, 0, 0, 0],
                'Sensitive to Drug': [0, 3, 4, 5, 2, 3],
                'Intermediate': [0, 0, 0, 0, 0, 0],
                'Resistant to Drug': [0, 3, 2, 1, 3, 3],
            })

    def test_paper_output(self):
        """Test paper sliced drug output"""
        drugs = self.assertJson('maps:data.drugs', filters=('paper',), data={'paper[]': 1})
        self.assertGraph(
            drugs,
            [None, 'BUMP', 'H2O', 'MEM', 'PIN', 'WAVE'], {
                'Unknown': [1, 0, 0, 0, 0, 0],
                'Sensitive to Drug': [0, 2, 4, 3, 5, 2],
                'Intermediate': [0, 0, 0, 0, 0, 0],
                'Resistant to Drug': [0, 7, 5, 6, 3, 7],
            })

    def test_map_output(self):
        """Test map sliced drug output"""
        drugs = self.assertJson('maps:data.drugs', filters=('map',), data={'map[]': ['DE', 'FR']})
        self.assertGraph(
            drugs,
            [None, 'BUMP', 'H2O', 'MEM', 'PIN', 'WAVE'], {
                'Unknown': [1, 0, 0, 0, 0, 0],
                'Sensitive to Drug': [0, 4, 4, 3, 3, 3],
                'Intermediate': [0, 0, 0, 0, 0, 0],
                'Resistant to Drug': [0, 4, 4, 4, 4, 4],
            })


class LineageData(BaseCase):
    """Test lineage data output."""
    def test_all_output(self):
        """Test not sliced lineage output"""
        lineages = self.assertJson('maps:data.lineages', field='children')
        self.assertEqual(lineages, [
            {'name': 'LA', 'color': 'rgb(48,129,189)', 'children': [], 'size': 4},
            {'name': 'LB', 'color': 'rgb(48,129,189)', 'children': [
                {'name': 'LB.1', 'color': 'rgba(48,129,189,0.8)', 'children': [], 'size': 4}
            ], 'size': 6},
            {'name': 'LC', 'color': 'rgb(48,129,189)', 'children': [], 'size': 6}
        ])

    def test_source_output(self):
        """Test source sliced lineage output"""
        source = ImportSource.objects.get(name='Import X')
        lineages = self.assertJson('maps:data.lineages', field='children',
                                   filters=('source',), data={'source[]': source.pk})
        self.assertEqual(lineages, [
            {'children': [], 'color': 'rgb(48,129,189)', 'name': 'LA', 'size': 1},
            {'children': [
                {'children': [], 'color': 'rgba(48,129,189,0.8)', 'name': 'LB.1', 'size': 1}
            ], 'color': 'rgb(48,129,189)', 'name': 'LB', 'size': 2},
            {'children': [], 'color': 'rgb(48,129,189)', 'name': 'LC', 'size': 2},
        ])

    def test_paper_output(self):
        """Test paper sliced lineage output"""
        lineages = self.assertJson('maps:data.lineages', field='children',
                                   filters=('paper',), data={'paper[]': 1})
        self.assertEqual(lineages, [
            {'children': [], 'color': 'rgb(48,129,189)', 'name': 'LA', 'size': 3},
            {'children': [
                {'children': [], 'color': 'rgba(48,129,189,0.8)', 'name': 'LB.1', 'size': 2}
            ], 'color': 'rgb(48,129,189)', 'name': 'LB', 'size': 2},
            {'children': [], 'color': 'rgb(48,129,189)', 'name': 'LC', 'size': 3},
        ])

    def test_map_output(self):
        """Test map sliced lineage output"""
        lineages = self.assertJson('maps:data.lineages', field='children',
                                   filters=('map',), data={'map[]': ['DE', 'RU']})
        self.assertEqual(lineages, [
            {'children': [], 'color': 'rgb(48,129,189)', 'name': 'LA', 'size': 4},
            {'children': [
                {'children': [], 'color': 'rgba(48,129,189,0.8)', 'name': 'LB.1', 'size': 4}
            ], 'color': 'rgb(48,129,189)', 'name': 'LB', 'size': 4},
            {'children': [], 'color': 'rgb(48,129,189)', 'name': 'LC', 'size': 5},
        ])

    def test_drug_output(self):
        """Test drug sliced lineage output"""
        lineages = self.assertJson('maps:data.lineages', field='children',
                                   filters=('drug',), data={'drug[]': ['WAVE', 'BUMP', 'PIN']})
        self.assertEqual(lineages, [
            {'children': [], 'color': 'rgb(48,129,189)', 'name': 'LA', 'size': 3},
            {'children': [
                {'children': [], 'color': 'rgba(48,129,189,0.8)', 'name': 'LB.1', 'size': 4}
            ], 'color': 'rgb(48,129,189)', 'name': 'LB', 'size': 6},
            {'children': [], 'color': 'rgb(48,129,189)', 'name': 'LC', 'size': 6},
        ])

class LocusListData(BaseCase):
    """Test list of locus names."""
    def test_all_output(self):
        """Test locus list is unsliced"""
        locus = GeneLocus.objects.get(gene_symbol='WA8')
        val = self.assertDataTable(
            'maps:data.locuses', names=('pk', 'str', 'start'), order=2,
            data={'genelocus[]': [locus.pk]},
            filters=()
        )
        vals = [unit['str'] for unit in val]
        self.assertEqual(vals, ['WA8', 'W1', 'W2', 'W3', 'W4', 'W5'])

class MutationsData(BaseCase):
    """Test mutations output, should only show mutations with resistances."""
    def test_all_output(self):
        """Test mutation list"""
        val = self.assertDataTable(
            'maps:data.mutations', order=0,
            names=('name', 'mode', 'gene_locus', 'nucleotide_position', 'strain_count'),
            data={'mutation[]': 'Mutation_020'},
            filters=()
        )
        vals = [unit['name'] for unit in val]
        self.assertEqual([unit['name'] for unit in val],\
            ['Mutation_020', 'Mutation_001', 'Mutation_002',
             'Mutation_003', 'Mutation_004', 'Mutation_005'])
        self.assertEqual([int(unit['strain_count']) for unit in val], [5, 10, 7, 3, 6, 11])


    def test_source_output(self):
        """Test mutations sliced by source"""
        sources = ImportSource.objects.filter(name__in=['Import Z', 'Import Y'])

        val = self.assertDataTable(
            'maps:data.mutations', order=0,
            names=('name', 'mode', 'gene_locus', 'nucleotide_position', 'strain_count'),
            data={'source[]': sources.values_list('pk', flat=True)},
            filters=('source',)
        )
        self.assertEqual([unit['name'] for unit in val],\
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_004', 'Mutation_005'])
        self.assertEqual([int(unit['strain_count']) for unit in val], [8, 4, 2, 5, 6])

    def test_paper_output(self):
        """Test mutations sliced by paper"""
        val = self.assertDataTable(
            'maps:data.mutations', order=0,
            names=('name', 'mode', 'gene_locus', 'nucleotide_position', 'strain_count'),
            data={'paper[]': 1},
            filters=('paper',)
        )
        self.assertEqual([unit['name'] for unit in val],\
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_004', 'Mutation_005'])
        self.assertEqual([int(unit['strain_count']) for unit in val], [5, 3, 1, 2, 7])

    def test_map_output(self):
        """Test mutations sliced by map"""
        val = self.assertDataTable(
            'maps:data.mutations', order=0,
            names=('name', 'mode', 'gene_locus', 'nucleotide_position', 'strain_count'),
            data={'map[]': ['DE', 'FR']},
            filters=('map',)
        )
        self.assertEqual([unit['name'] for unit in val],\
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_005', 'Mutation_006'])
        self.assertEqual([int(unit['strain_count']) for unit in val], [4, 2, 1, 5, 3])

    def test_locus_output(self):
        """Test mutations sliced by locus"""
        loci = GeneLocus.objects.filter(gene_symbol__in=['W9', 'WA8', 'W4'])
        val = self.assertDataTable(
            'maps:data.mutations', order=0,
            names=('name', 'mode', 'gene_locus', 'nucleotide_position', 'strain_count'),
            data={'genelocus[]': loci.values_list('pk', flat=True)},
            filters=('genelocus',)
        )
        self.assertEqual([unit['name'] for unit in val],\
            ['Mutation_006', 'Mutation_007', 'Mutation_016', 'Mutation_017', 'Mutation_018'])
        self.assertEqual([int(unit['strain_count']) for unit in val], [6, 6, 4, 5, 1])

    def test_drug_output(self):
        """Test mutations sliced by drug"""
        val = self.assertDataTable(
            'maps:data.mutations', order=0,
            names=('name', 'mode', 'gene_locus', 'nucleotide_position', 'strain_count'),
            data={'drug[]': ['WAVE', 'PIN', 'MEM']},
            filters=('drug',)
        )
        self.assertEqual([unit['name'] for unit in val],\
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_004', 'Mutation_005'])
        self.assertEqual([int(unit['strain_count']) for unit in val], [9, 7, 3, 6, 10])


class MutationResistanceData(BaseCase):
    """Test mutation resistance data."""
    def test_all_output(self):
        """Test mutation list"""
        mutations = self.assertJson('maps:data.mutation', filters=('mutation',),\
            data={'mutation[]': ['Mutation_020', 'Mutation_001', 'Mutation_002', 'Mutation_003']})
        self.assertGraphTotals(
            mutations,
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_020'], {
                # KeyName:  (Numerators, Denominators),
                # CONFIRMED DATA
                'Sensitive': ([1, 0, 1, 1], [2, 2, 2, 2]),
                'Other Drug Resistant': ([0, 0, 0, 0], [-1, -1, -1, -1]),
                'Multi Drug Resistant': ([7, 6, 2, 2], [13, 13, 13, 13]),
                'Extensively Drug Resistant': ([2, 1, 0, 2], [5, 5, 5, 5]),
            })

    def test_source_output(self):
        """Test mutations sliced by source"""
        sources = ImportSource.objects.filter(name__in=['Import Z', 'Import Y'])

        mutations = self.assertJson('maps:data.mutation', filters=('source', 'mutation'), data={
            'source[]': sources.values_list('pk', flat=True),
            'mutation[]': ['Mutation_020', 'Mutation_001', 'Mutation_002', 'Mutation_003'],
        })
        self.assertGraphTotals(
            mutations,
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_020'], {
                # CONFIRMED DATA
                'Sensitive': ([1, 0, 1, 1], [2, 2, 2, 2]),
                'Other Drug Resistant': ([0, 0, 0, 0], [-1, -1, -1, -1]),
                'Multi Drug Resistant': ([6, 4, 1, 1], [10, 10, 10, 10]),
                'Extensively Drug Resistant': ([1, 0, 0, 2], [2, 2, 2, 2]),
            })

    def test_paper_output(self):
        """Test mutations sliced by paper"""
        mutations = self.assertJson('maps:data.mutation', filters=('paper', 'mutation'), data={
            'paper[]': 1,
            'mutation[]': ['Mutation_020', 'Mutation_001', 'Mutation_002', 'Mutation_003'],
        })
        self.assertGraphTotals(
            mutations,
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_020'], {
                # CONFIRMED DATA
                'Sensitive': ([0, 0, 0, 0], [-1, -1, -1, -1]),
                'Other Drug Resistant': ([0, 0, 0, 0], [-1, -1, -1, -1]),
                'Multi Drug Resistant': ([3, 2, 1, 1], [5, 5, 5, 5]),
                'Extensively Drug Resistant': ([2, 1, 0, 2], [5, 5, 5, 5]),
            })

    def test_map_output(self):
        """Test mutations sliced by map"""
        mutations = self.assertJson('maps:data.mutation', filters=('map', 'mutation'), data={
            'map[]': ['DE', 'FR'],
            'mutation[]': ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_005'],
        })
        self.assertGraphTotals(
            mutations,
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_005'], {
                # CONFIRMED DATA
                'Sensitive': ([0, 0, 0, 0], [1, 1, 1, 1]),
                'Other Drug Resistant': ([0, 0, 0, 0], [-1, -1, -1, -1]),
                'Multi Drug Resistant': ([3, 2, 1, 4], [7, 7, 7, 7]),
                'Extensively Drug Resistant': ([1, 0, 0, 1], [1, 1, 1, 1]),
            })

    def test_drug_output(self):
        """Test mutations sliced by drug"""
        mutations = self.assertJson('maps:data.mutation', filters=('drug', 'mutation'), data={
            'drug[]': ['PIN'],
            'mutation[]': ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_004'],
        })
        self.assertGraphTotals(
            mutations,
            ['Mutation_001', 'Mutation_002', 'Mutation_003', 'Mutation_004'], {
                # CONFIRMED DATA
                # s: BOB2,BOB6,BOB9,BAB2,BAB3,BAB4,BAB7,BAB8
                'Sensitive to Drug': ([6, 3, 1, 2], [8, 8, 8, 8]),
                'Intermediate': ([0, 0, 0, 0], [-1, -1, -1, -1]),
                # r: BOB3,BOB4,BOB5,BOB8,BAB1,BAB5,BAB6,BBB1
                'Resistant to Drug': ([2, 3, 1, 3], [8, 8, 8, 8]),
            })

    def test_multi_drug_output(self):
        """Test multiple drugs selected"""
        mutations = self.assertJson('maps:data.mutation', filters=('drug', 'mutation'), data={
            'drug[]': ['PIN', 'WAVE'],
            'mutation[]': ['Mutation_001', 'Mutation_002',],
        })
        print(mutations)
        self.assertGraphTotals(
            mutations,
            ['Mutation_001 (PIN)', 'Mutation_001 (WAVE)', 'Mutation_002 (PIN)', 'Mutation_002 (WAVE)'], {
                # PIN/r: BOB3,BOB4,BOB5,BOB8,BAB1,BAB5,BAB6,BBB1
                # WAVE/r: BOB1,BOB2,BOB6,BOB8,BAB2,BAB3,BAB4,BAB5,BAB8,BBB1,BBB2
                'Sensitive to Drug': ([6, 3, 3, 4], [8, 7, 8, 7]),
                'Intermediate': ([0, 0, 0, 0], [-1, -1, -1, -1]),
                # PIN/s: BOB2,BOB6,BOB9,BAB2,BAB3,BAB4,BAB7,BAB8
                # WAVE/s: BOB3,BOB4,BOB5,BOB9,BAB1,BAB6,BAB7
                'Resistant to Drug': ([2, 6, 3, 3], [8, 11, 8, 11]),
            })
