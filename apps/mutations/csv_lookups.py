#
# Copyright (C) 2017-2018  Dr. Maha Farhat
#               2017-2018  Martin Owens (Author)
#                    2018  MIT, TESS
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
Opens CSV/TSV files and indexes their contents into tables which can then be
used to look data up.

Features data conversion, sub table support and column naming.
"""

import sys
import logging

from collections import OrderedDict
from csv import reader as CsvReader

class BaseLookup(dict):
    """
    Master lookup structure opens TSV files and indexes columns against
    each other. Multiple files can be concatinated together and column names
    can be replaced to make things consistant.

    dat = Summary('filename.tsv', column_wrong='column_correct')
    dat.load_file('second_filename.tsv')

    """
    @property
    def key(self):
        raise NotImplementedError("Key lookup column must be specified in parent class")

    null='\\N'
    delimiter='\t'
    reverse = None
    defaults = {}
    replaces = {}
    types = {}
    sub_tables = None

    def __init__(self, filename=None, **kw):
        if filename is not None:
            self.load_file(filename, **kw)

    def load_file(self, filename, **syn):
        """Loads a file into this summary"""
        with open(filename, 'r') as fhl:
            # Open the TSV file and get the header as the first row.
            rows = list(CsvReader(fhl, delimiter=self.delimiter))
            header = rows.pop(0)

            # Rename the header to more useful names
            for a, b in syn.items():
                if a in header and b not in header:
                    header[header.index(a)] = b

            # Stow header for usefulness
            self.header = header

            # Add each of the rows to the internal lookups
            for line, row in enumerate(rows):
                if not line:
                    # Ignore empty rows
                    continue
                try:
                    self.append(dict(zip(header, row)))
                except ValueError as err:
                    logging.warning("%s:%d %s" % (filename, line+2, str(err)))

    def set_type(self, key, value):
        """Some of the data in the file isn't a string, parse it"""
        return value if value is None else self.types.get(key, str)(value)

    def append(self, row):
        """Append one row from the summary file."""
        key = row.get(self.key, None)
        if key is None:
            logging.warning("Ignoring row with no key {}".format(self.key))

        if key not in self:
            self[key] = row
        else:
            return self.append_tables(row, key)

        for name in row:
            # Create the initial dictionary container for this lookup
            if not hasattr(self, name):
                setattr(self, name, OrderedDict())

            # Use default values for some rows if the column was empty
            if row[name] == self.null:
                row[name] = None
            if not row[name]:
                row[name] = self.defaults.get(name, row[name])

            # Set the type if needed, useful for later use
            row[name] = self.set_type(name, row[name])

            # And then mutate the value with replacement characters if needed
            if isinstance(row[name], str):
                for a, b in self.replaces.get(name, ()):
                    row[name] = row[name].replace(a, b)

            # The primary key often has a reverse lookup to a second column
            lookup = self.reverse if name == self.key else self.key

            # Skipped columns are ignored when requested
            if lookup is None or name is None:
                continue

            key = self.set_type(lookup, row[lookup])

            # Load the internal dictionary read for use
            loc = getattr(self, name)

            # Protect the dictionaries from duplicate key errors
            if key in loc:
                raise ValueError("Duplicate reverse index %s=%s = %s|%s" % (
                    lookup, key, row[name], loc[key]))

            # Save the data into the dictionary
            loc[key] = row[name]

    def append_tables(self, row, key):
        """Pop out any sub tables as needed into their sub-structures"""
        oldrow = self[key]
        for sdef in (self.sub_tables or []):
            if callable(sdef):
                sdef = sdef(self.header, row, key)

            if isinstance(sdef, str): # Creates a single list
                if not isinstance(oldrow[sdef], list):
                    oldrow[sdef] = [oldrow[sdef]]
                oldrow[sdef].append(row.pop(sdef))

            elif isinstance(sdef, (list, tuple)): # Creates a table
                def _parse(target):
                    subrow = [target.pop(n) for n in sdef]
                    if len(subrow) > 2:
                        return {subrow[0]: dict(zip(sdef[1:], subrow[1:]))}
                    return {subrow[0]: subrow[1]}

                if not isinstance(oldrow[sdef[0]], dict):
                    oldrow[sdef[0]] = _parse(oldrow)
                oldrow[sdef[0]].update(_parse(row))

        # Check the rest of the values for sanity, allow null values in rows.
        for name in row:
            if row[name] is None:
                continue
            elif self[key][name] is None:
                self[key][name] = row[name]
            elif row[name] != self[key][name]:
                raise KeyError("Duplicate key with different data: %s=%s %s=(%s!=%s)" % (self.key, key, name, self[key][name], row[name]))


class Lookup(BaseLookup):
    """A basic lookup"""
    auto_delim = {'tsv': '\t', 'csv': ','}

    def __init__(self, filename, delimiter=None, key=None, sub_tables=None, **kw):
        self._key = key

        if sub_tables is not None:
            self.sub_tables = sub_tables

        if delimiter is None:
            ext = filename.rsplit('.', 1)[-1]
            delimiter = self.auto_delim.get(ext, None)

        if delimiter is not None:
            self.delimiter = delimiter

        super(Lookup, self).__init__(filename, **kw)

    @property
    def key(self):
        if self._key:
            return self._key
        if self.header:
            return self.header[0]
        return super(Lookup, self).key

