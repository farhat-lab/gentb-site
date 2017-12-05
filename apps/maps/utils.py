#
# Copyright (C) 2017  Dr. Maha Farhat
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
Data management for maps, basic functions.
"""

from collections import defaultdict, OrderedDict

def to(method):
    """Turn generators into objects, method can be a type, obj or function"""
    def _outer(f):
        def _inner(*args, **kw):
            return method(f(*args, **kw))
        return _inner
    return _outer

class OrderlyDict(OrderedDict):
    """Special OrderedDict for pre-processing"""
    def __init__(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            args = list(args)
            if args[0] and not isinstance(args[0][0], (list, tuple)):
                args[0] = zip(args[0], args[0])
        super(OrderlyDict, self).__init__(*args, **kwargs)

class GraphData(defaultdict):
    """Format three columns into a format suitable for d3 graphs"""
    def __init__(self, qs, x, y, z, trim=False):
        super(GraphData, self).__init__(lambda: defaultdict(int))
        self.keys = defaultdict(OrderedDict)
        self.trims = {'x': trim, 'y': trim, 'z': trim}

        for dd in qs:
            # Collapse multiple fields into categories
            if isinstance(x, list):
                for tx in x:
                    self[tx][dd[tx]] += dd[y]
            # Or take categories from one field
            elif dd[y] > 0:
                self.keys['x'][dd[x]] = dd[x]
                self[dd.get(z, None)][dd[x]] += dd[y]

    def set_axis(self, axis, keys=None, trim=None):
        """
        Very important function for defining how the data is collated.

          axis - The axis we are setting, this can be
            x: often called name or series
            y: ofteh the value or magnatude
            z: usually the category, not always used.

          keys - A list of tuples containing an ordered match between
                 values found in the database and display values. Extra
                 keys can be added, even if they don't appear in the db.

          trim - Sets this axis to be trimmed, this is where the value is
                 empty or none and it is then exluded from output.
                   x - The name is None or empty string, so ignored, setting
                       the keys to ('Something', None) will exclude them.
                   y - The value is zero or None, so ignored.
                   z - The category contains no rows
                       (maybe because x and y trimming)
        """
        if trim is not None:
            self.trims[axis] = trim
        if keys is not None:
            self.keys[axis].update(OrderlyDict(keys))
        return self

    def get_z_cats(self):
        return self.keys.get('z', OrderlyDict(list(self)))

    def get_x_cols(self, cat):
        cols = OrderlyDict(list(self[cat]))
        ret = self.keys.get('x', cols).copy()
        ret.update(cols)
        return ret

    def get_y_value(self, cat, col):
        total = self.keys['y'].get(cat)
        value = self[cat][col]
        if total is not None:
            if total is 0:
                return 0, value, total
            else:
                return value / float(total), value, total
        return value, value, -1

    def get_values(self, cat):
        """Generator returning all values in this category"""
        for col, name in self.get_x_cols(cat).items():
            y, value, total = self.get_y_value(cat, col)
            if self.is_trimmed(name, value, cat):
                yield {
                    "y": y,
                    "x": name,
                    "col": col,
                    "value": value,
                    "total": total,
                }

    def is_trimmed(self, name, value, cat):
        x, y = self.trims['x'], self.trims['y']
        if isinstance(x, (list, tuple)):
            x = cat in x
        if isinstance(y, (list, tuple)):
            y = cat in y
        return (value or not y) and (name or not x)

    @to(list)
    def to_graph(self):
        """Make square structure and convert from defaultdict to OrderedDict"""
        for cat, name in self.get_z_cats().items():
            values = list(self.get_values(cat))
            if values or not self.trims['z']:
                yield {"key": name, "values": values}


