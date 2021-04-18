#
# Copyright (C) 2021 Maha Farhat
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
Test for mutations
"""

import json
from collections import Counter
from itertools import zip_longest

from django.test import TestCase
from extratest.base import ExtraTestCase

from .models import Drug, StrainMutationCache
from ..maps.models import Country

def set_pk(model, pk, **kw):
    con = model.objects.get(**kw)
    con.delete()
    con.pk = 3
    con.save()

class CountCacheTest(TestCase):
    """Test caching of count queries"""
    fixtures = ['test-genetics', 'test-maps', 'test-strains']
    maxDiff = 30000

    def setUp(self):
        set_pk(Country, 3, iso2='DZ')
        set_pk(Drug, 5, code='H2O')

    def assertCacheKey(self, key, **kwargs):
        """Assert the cache is correct"""
        expected = {}
        for x in ('importer', 'paper', 'country', 'lineage'):
            expected[x] = key if x in kwargs else None
        self.assertEqual(StrainMutationCache.key_from_inputs(**kwargs), expected)

    def test_caching_keys(self):
        """Can generate caching keys"""
        self.assertCacheKey('')
        self.assertCacheKey('AQ==', importer=[1])
        self.assertCacheKey('Ag==', paper=[2])
        self.assertCacheKey('AgM=', paper=[2, 3])
        self.assertCacheKey('AAM=', country=['DZ'])
        #self.assertCacheKey('Aw==', drug=['H2O'])
        self.assertCacheKey('LGI=', lineage=['B'])

    def test_matrix_call(self):
        self.ret = []
        self.count = 0
        def fun(**kwargs):
            self.ret.append([
                str(item) for item in kwargs.values()
            ])
            self.count += 1
            if self.count > 10:
                return False
        StrainMutationCache.matrix_call(fun)
        self.assertEqual(self.ret, [
            ['None', 'None', 'None', 'None'],
            ['None', 'None', 'None', 'A'],
            ['None', 'None', 'None', 'B'],
            ['None', 'None', 'None', 'B.1'],
            ['None', 'None', 'None', 'C'],
            ['None', 'None', 'Afghanistan', 'None'],
            ['None', 'None', 'Afghanistan', 'A'],
            ['None', 'None', 'Afghanistan', 'B'],
            ['None', 'None', 'Afghanistan', 'B.1'],
            ['None', 'None', 'Afghanistan', 'C'],
            ['None', 'None', 'Åland Islands', 'None'],
        ])

