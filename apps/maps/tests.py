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

import os
import time
import signal
import tempfile

from django.test import TestCase, override_settings
from autotest.base import ExtraTestCase

from .utils import OrderlyDict, OrderedDict, GraphData

class UtilsTest(TestCase):
    def test_orderly(self):
        a = OrderlyDict(OrderedDict([('a', 'a'), ('b', 'b'), ('c', 'c')]))
        self.assertEqual(a, OrderlyDict(['a', 'b','c']))

        b = OrderedDict([('a', 1), ('b', 2), ('c', 3)])
        self.assertEqual(b, OrderlyDict([('a', 1), ('b', 2), ('c', 3)]))

    def test_graph(self):
        raise NotImplementedError("Need to be written")
