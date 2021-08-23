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

COUNTRY_MAP = {
    'unknown': None,
    'UK': 'United Kingdom',
    'United Kingdom of Great Britain and Northern Ireland': 'United Kingdom',
    'United States of America': 'United States',
    'South Afica': 'South Africa',
    'Ivory Coast': 'CIV',
    "Côte d'Ivoire": 'CIV',
    'Netherland': 'Netherlands',
    'Camerun': 'Cameroon',
    'Trinidad & Tobago': 'Trinidad and Tobago',
    'Kazakstan': 'Kazakhstan',
    'Kazachstan': 'Kazakhstan',
    'Brasil': 'Brazil',
    'Azerbaidjan': 'Azerbaijan',
    'Marocco': 'Morocco',
    'DR Congo': 'Democratic Republic of the Congo',
    'RD Congo': 'Democratic Republic of the Congo',
    'Guinea Eq.': 'Equatorial Guinea',
    'Myanmar': 'Burma',
    'Netherlands Antilles': 'Sint Maarten',
    'Guinea-Conakry': 'Guinea',
    'Czechia': 'Czech Republic',
    'China /Tibet': 'Tibet',
    'Carribean': 'Aruba',
    'Russian Federation': 'Russia',
    'Philipines': 'Philippines',
    'Comoro Islands': 'Comoros',
    'South Korea N': 'Korea, Republic of',
    'Korea': 'Korea, Republic of',
    'Republic of Korea': 'Korea, Republic of',
    'Eswatini': 'Swaziland',
    'China, Hong Kong SAR': 'Hong Kong',
}

CITY_MAP = {
    'Karakalpakstan': 'Nukus', # Region
    'New York City': 'New York',
    'KwaZulu-Natal': 'Durban',
    'Kwazulu Natal': 'Durban',
    'South Carolina': 'Columbia',

    'isolated in SF Lineage3A': None,
    'CDC': None,
}

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
    def __init__(self, qset, x_axis, y_axis, z_axis, trim=False, filter_label=None):
        super(GraphData, self).__init__(lambda: defaultdict(int))
        self._filter_label = filter_label
        self.keys = defaultdict(OrderedDict)
        self.trims = {'x': trim, 'y': trim, 'z': trim}

        for row in qset:
            # Collapse multiple fields into categories
            if isinstance(x_axis, list):
                for col_name in x_axis:
                    x_col = self.filter_label('x', row[col_name])
                    y_col = self.filter_label('y', row[y_axis])
                    z_col = self.filter_label('z', col_name)
                    self[z_col][x_col] += y_col
            # Or take categories from one field
            elif row[y_axis] > 0:
                x_col = self.filter_label('x', row[x_axis])
                y_col = self.filter_label('y', row[y_axis])
                z_col = self.filter_label('z', row.get(z_axis, None))
                self.keys['x'][x_col] = x_col
                self[z_col][x_col] += y_col

    def filter_label(self, axis, label):
        """Apply any filters to the axis labels"""
        if label is None or not self._filter_label:
            return label
        return self._filter_label(axis, label)

    def set_axis(self, axis, keys=None, trim=None):
        """
        Very important function for defining how the data is collated.

          axis - The axis we are setting, this can be
            x: often called name or series
            y: often the value or magnatude
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
        # I believe the best output is to use the names from keys
        # and the orders from cols. But the is a trade-off
        cols = OrderlyDict(list(self[cat]))
        cols.update(self.keys.get('x', cols))
        return cols

    def get_y_value(self, cat, col):
        total = self.keys['y'].get(cat)
        value = self[cat][col]
        if total is not None:
            if total is 0:
                return 0, value, total
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

class Sdict(OrderedDict):
    """
    Decodes query dictionaries (from django's QueryDict) into multi-level python structures.

    For example, order__0__foo = 'desc' becomes {'order': [{'foo': 'desc'}]}
    """
    _filter_key = staticmethod(lambda key: key)
    def __init__(self, querydict=None):
        super(Sdict, self).__init__()
        if querydict is not None:
            self.update(querydict)

    def update(self, querydict):
        """Add the dictionary, parsing out multi-level values"""
        for (key, value) in querydict.items():
            self[self._filter_key(key)] = value
        self.toarrays()

    def toarrays(self):
        """Looks at the keys in this dictionary and returns it as an array
        if all the keys are sequential integers."""
        for key, value in self.items():
            if isinstance(value, Sdict):
                self[key] = self[key].toarrays()
        try:
            ints = sorted([(int(key), key) for key in self])
            if not ints:
                return []
            if tuple(next(iter(zip(*ints)))) != tuple(range(len(ints))):
                raise ValueError("Integers aren't sequential.")
            return [self[key] for _, key in ints]
        except ValueError:
            return self

    def __setitem__(self, key, value):
        if isinstance(key, str):
            keys = key.split('__')
        elif isinstance(key, (list, tuple)):
            keys = key
        else:
            raise TypeError("Unknown key type: {}".format(type(key).__name__))

        if not keys:
            return None # Don't set empty keys

        key = keys[0]

        if len(keys) == 1:
            return super(Sdict, self).__setitem__(key, value)

        if key not in self or not isinstance(self[key], OrderedDict):
            super(Sdict, self).__setitem__(key, type(self)())

        self[key][keys[1:]] = value
        return value

class Jdict(Sdict):
    """
    Decodes data tables ajax request information into python. (see Sdict for details)

    For example, order[0][foo] = 'desc' becomes {'order': [{'foo': 'desc'}]}

    All keys and values are treated genertically with no foreknowlege about dataTables.
    """
    @staticmethod
    def _filter_key(key):
        return key.replace('][', '__').replace('[', '__').replace(']', '')

def many_lookup(model, local_filter, id_field, ret_field='pk__in'):
    """
    When looking up many-to-many values, we need to not duplicate rows
    in the output when it comes to counts. This provides us with a way
    of joining things up.
    """
    def _inner(values):
        return ret_field,\
            model.objects.filter(**{local_filter+'__in': values})\
                         .values_list(id_field, flat=True)
    return _inner

def adjust_coords(geom):
    """
    Attempt to adjust the generate coordinates to fit leaflet map projection.
    """
    for item in geom['coordinates']:
        _adjust(item)
    return geom

def _adjust(lst):
    if not isinstance(lst, list):
        raise ValueError("Can't adjust non-list")
    if lst and isinstance(lst[0], list):
        for item in lst:
            _adjust(item)
    elif len(lst) == 2 and isinstance(lst[0], (float, int)):
        lst[0], lst[1] = lst[1], lst[0]

