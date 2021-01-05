#
# Copyright (C) 2018-2019 - Dr. Maha Farhat
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
Match and extract information from snp names and other encoded information.
"""

import os
import re
import csv
import json
import sys

from functools import reduce
from collections import defaultdict, OrderedDict
from operator import or_
from datetime import date

from django.db.models import Q

MONTHS = ['', 'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
DATE_FORMATS = [
    r'^(?P<year>\d{4})-(?P<mm>\d{2})-(?P<dd>\d{2})$',
    r'^(?P<year>\d{2,4})$',
    r'^(?P<mm>\d{1,2})/(?P<year>\d{2,4})$',
    r'^(?P<mon>[A-Za-z]{3})/(?P<year>\d{2,4})$',
    r'^(?P<mm>\d{1,2})/(?P<dd>\d{1,2})/(?P<year>\d{2,4})$',
    r'^(?P<dd>\d{2})(?P<mon>[A-Za-z]{3})(?P<year>\d{2,4})$',
]

def re_match_raw(re_list, string):
    """
    Match the first regular expression in a list of regular expressions
    and return the groupdict of the matching expression.

    Will replace your list items with re objects in place! for caching.
    """
    string = str(string).strip()
    for (x, rex) in enumerate(re_list):
        # Compile only once per run
        if isinstance(rex, str):
            rex = re.compile(rex)
            re_list[x] = rex
        match = rex.search(string)
        if match is not None:
            return x, match
    raise ValueError("Couldn't match string: %s" % str(string))

def re_match(re_list, string):
    """Returns the groupdict for re_match_raw"""
    return re_match_raw(re_list, string)[-1].groupdict()

def re_match_dict(re_dict, string):
    """
    Like re_match(), but instead of a list, it takes a dictionary. The first
    matching regular expression then has it's value populated with the
    match's group dictionary as a template. So for example:

    re_match_dict({r'P(\\d+)': 'A%s'}, 'P45') == 'A45'

    It will also take a list of tuples for ordered matching.
    """
    string = str(string).strip()
    re_dict = re_dict.items() if isinstance(re_dict, dict) else re_dict
    for (rex, value) in re_dict:
        # Compile only once per run
        if isinstance(rex, str):
            rex = re.compile(rex)
        match = rex.search(string)
        if match is not None:
            if '%' in value:
                return value % match.groupdict()
            return value
    raise ValueError("Couldn't match string: %s" % str(string))

def pop_all(data, *args):
    """Remove all keys from a dictionary"""
    for key in args:
        data.pop(key, None)


def get_date(string):
    """Attempt to parse an unknown and variable formated date"""
    if string in [None, '', u'', '.', 0, False]:
        return None
    (year, month, day) = (1970, 1, 1)
    bits = re_match(DATE_FORMATS, string)
    if 'year' in bits:
        year = bits['year']
        if len(year) == 2:
            if int(year) > 50:
                year = int(year) + 1900
            else:
                year = int(year) + 2000
        elif len(year) == 4:
            year = int(year)
        else:
            raise ValueError("Invalid year: %s" % str(bits))
    if 'mon' in bits:
        month = MONTHS.index(bits['mon'].strip('0').lower())
    elif 'mm' in bits:
        month = int(bits['mm'])
    if 'dd' in bits:
        day = int(bits['dd'])
    return date(year, month, day)

CODING = r'(?P<coding>(?P<cref>[ACTG])(?P<cpos>\d+)(?P<cver>[ACTG]))'
CODINGS = r'(?P<coding>(?P<cref>[ACTG]+)(?P<cpos>\d+)-(?P<cpos2>\d+)(?P<cver>[ACTG]+))'
NONCODE = r'(?P<noncode>promoter|inter)'
AMINO = r'(?P<amino>(?P<aref>[A-Z\*])(?P<apos>\d+)(?P<aver>[A-Z\*]))'
AMINOS = r'(?P<amino>(?P<aref>[A-Z\*]+)(?P<apos>\d+)-(?P<apos2>\d+)(?P<aver>[A-Z\*]+))'

SNP = r'^SNP_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_' + CODING
LSP = r'^LSP_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_' + CODINGS
GENE = r'(?P<gene>[a-zA-Z\d\-_\.]+)|(?P<rgene>rr[sl]))\'?$'

MUTATION_RE = [
    SNP + r'_((' + AMINO + r'|' + NONCODE + ')[-_]' + GENE,
    LSP + r'_((' + AMINOS + r'|' + NONCODE + r')[-_]' + GENE,
    r'^(?P<mode>(INS|DEL))_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_(i|d|\.|i\.)?'\
        r'(?P<codes>(?P<cpos>[\d\-]+)(?P<cref>[ATGC]*))_((?P<noncode>promoter|inter|\d+)_)?'\
        r'(?P<gene>[a-zA-Z\d\-_\.]+?)(_' + AMINO + r')?\'?$',
    # These are older SNP names and should probably be converted
    SNP + r'_(?P<gene>(Rv\w+|rrl|rrs))',
    SNP + r'_PE_(?P<gene>[a-zA-Z\d\-]+)_' + AMINO + r'$',
    SNP + r'_(?P<gene>[a-zA-Z\d\-]+)_' + AMINO + r'$',
    SNP + r'_?[\.A-Z](\d+)[\.A-Z]_(?P<gene>[a-zA-Z\d\-_\.]+)$',
    SNP + r'_(?P<gene>[a-zA-Z\d\-_]+)$',
]

def match_snp_name(name):
    """
    Tries to regex match the snp name and returns a dictionary of info decoded from the SNP name.
    """
    return re_match(MUTATION_RE, name)

def match_snp_half(name):
    """
    Backup attempt to match the ntpos only (for cross reference of bad names)
    """
    return re_match([SNP], name)

def match_snp_name_raw(name):
    """Same as above but raw values returned"""
    return re_match_raw(MUTATION_RE, name)

SYN_A = {
    'I': 'Integenic',
    'P': 'Promoter',
    'C': 'Coding',
    'N': 'Noncoding',
}
SYN_B = {
    'D': 'Deletion',
    'I': 'Insertion',
    'F': 'Frame shift',
    'S': 'Silent or Synonymous',
    'Z': 'Stop codon, nonsense mutation',
}

def info_mutation_format(mutation):
    """
    Processes a mutation name into a dictionary of parsing statements
    """
    # Process the snp name into it's parts
    index, match = match_snp_name_raw(mutation)
    items = [(name,) + match.span(name) for name in match.groupdict()]

    # Sort by the start of the match
    items.sort(key=lambda i: i[1])

    # From the end to the start, hilight the text.
    ret = []
    snp = OrderedDict()
    last_end = 0
    for (name, start, end) in items:
        if start == end:
            continue
        if start > last_end:
            ret.append(mutation[last_end:start])
        value = mutation[start:end]
        snp[name] = value
        ret.append('<span class="match %s">%s</span>' % (name, value))
        last_end = end

    # Prepend any remaining to the highlighted version.
    if last_end != len(mutation) - 1:
        ret.append(mutation[last_end:])

    return ''.join(ret), MUTATION_RE[index].pattern, snp


def unpack_mutation_format(name):
    """
    Processes a mutation into it's three main parts:

     * Index - Optional order of the mutation
     * Gene - The gene involved in this mutation
     * Mutation - Returned Mutation, usually just cleaned.

    """
    index = None
    if " " in name:
        index, name = name.split(" ", 1)
        try:
            index = int(index)
        except:
            raise ValueError("Optional sort index should be a number.")
    if index is None:
        index = 0

    snp = match_snp_name(name)
    orig = snp.get('gene', None)
    if orig is None:
        orig = snp.get('rgene', None)
    if orig is None:
        raise ValueError("Can't find gene name in %s" % name)
    # Normalise the inter-gene seperator
    gene = orig.replace('.', '-').replace('_', '-').strip("'").strip("-")

    # Put the inter-gene back into the mutation name
    name = name.replace(orig, gene)
    if index > 1:
        # These are older formats and should be re-formatted
        name = generate_mutation_name(**snp)

    if snp['syn'] == 'P':
        if snp.get('noncode', '') != 'promoter':
            _ = ValueError("Promoter doesn't specify 'promoter' part in %s" % name)
        return (index, 'promoter ' + gene, name)
    if snp['syn'] == 'I':
        if snp.get('noncode', '') != 'inter':
            _ = ValueError("Integenic doesn't specify 'inter' part")
        return (index, 'intergenic ' + gene, name)
    if snp['syn'] in ['CN', 'CD', 'CF', 'CI', 'CS', 'CZ', 'N', 'ND', 'NI', 'NF']:
        return (index, gene, name)
    raise ValueError("Must be promoter, intergenic or CN, CD, CF, CI, CS, CZ or N, ND, NI, NF")


def generate_mutation_name(mode='SNP', **snp):
    """Takes the properties and attempts to generate a name"""
    snp['mode'] = mode
    snp['gene'] = snp.get('gene', snp.pop('rgene', None))
    snp['coding'] = snp.get('coding', snp.pop('codes', None))
    if 'amino' in snp:
        if 'noncode' in snp:
            return '%(mode)s_%(syn)s_%(ntpos)s_%(coding)s_%(noncode)s_%(gene)s_%(amino)s' % snp
        return '%(mode)s_%(syn)s_%(ntpos)s_%(coding)s_%(amino)s_%(gene)s' % snp
    return '%(mode)s_%(syn)s_%(ntpos)s_%(coding)s_%(noncode)s_%(gene)s' % snp


class defaultlist(defaultdict): # pylint: disable=invalid-name
    """Like a defaultdict, but generates a list of items with the same key"""
    def __init__(self, generator):
        super(defaultlist, self).__init__(list)
        for key, value in generator:
            self[key].append(value)

    def flatten(self, method, **kw):
        """Trun each list value into a single value based on method"""
        for key, value in self.items():
            self[key] = method(value, **kw)

    def re_key(self, new_key):
        """Assuming each value is a dictionary, re-key the dictionary"""
        for key in list(self):
            value = self.pop(key)
            self[value.pop(new_key)].append(value)

class FileNotFound(IOError):
    """File not found error"""

class FieldsNotFound(KeyError):
    """Fields not found error"""

def csv_merge(fhl, **kw):
    """Merge csv header into each row for dictionary output"""
    reader = csv.reader(fhl, **kw)
    header = next(reader)
    yield header
    for row in reader:
        yield dict(zip(header, row))

LOADERS = {
    'json': lambda fhl: (json.loads(fhl.read()), False),
    'csv': lambda fhl: (csv_merge(fhl, delimiter=','), True),
    'tsv': lambda fhl: (csv_merge(fhl, delimiter='\t'), True),
}

def file_generator(*required, **_):
    """Decorate a method that uses csv data"""
    required = set(required)
    def _check(header, filename):
        missing = (header ^ required) & required
        if missing:
            raise FieldsNotFound("Fields '%s' missing from file %s" %\
                ("', '".join(missing), filename))

    def _outer(func):
        def _inner(*args, **kw):
            #args = list(args)
            # Handle cases of 'self' being the first argument
            index = int(not isinstance(args[0], str) and len(args) > 1)
            filename = args[index]

            # Make sure the file really exists.
            if not os.path.isfile(filename):
                raise FileNotFound("File '%s' Not Found" % filename)

            # Get the right content unpacker
            _loader = LOADERS.get(kw.get('loader', None))\
                   or LOADERS.get(filename.rsplit('.', 1)[-1])
            if _loader is None:
                raise TypeError("Can't parse '%s' unknown type." % filename)

            with open(filename, 'r') as fhl:
                (rows, header) = _loader(fhl)
                if 'status' in kw:
                    rows = StatusBar(kw.pop('status'), len(rows), rows, True)
                if header:
                    _check(set(next(rows)), filename)

                for row in rows:
                    value = func(*(args[:index] + (row,) + args[index+1:]), **kw)
                    if not header:
                        _check(set(row), filename)
                    if value is not None:
                        yield value
        return _inner
    return _outer

def json_generator(func):
    """Decorate a method that processes one row in a json filename list"""
    # Backwards compatible
    return file_generator(loader='json')(func)

def to(method): # pylint: disable=invalid-name
    """Turn generators into objects, method can be a type, object or function"""
    def _outer(func):
        def _inner(*args, **kw):
            return method(func(*args, **kw))
        return _inner
    return _outer

@to(dict)
@json_generator
def json_dict(row, key_id):
    """Creates a dictionary using the key_id, removes key_id from row"""
    return (row.pop(key_id), row)

@to(defaultlist)
@json_generator
def json_dictlist(row, key):
    """Creates a list of dictionaries using the key as a template"""
    return (key % row, row)

def get_bool(string):
    """Forces the datum to a boolean"""
    return str(string).upper() not in ['N', 'NO', 'FALSE', 'F', '0', 'NONE']

def get_int(datum):
    """Forces the datum into an integer"""
    if isinstance(datum, str):
        try:
            return int(datum)
        except ValueError:
            if len(datum) == 1:
                return ord(datum)
    return datum

TR_METHODS = {
    'date': get_date,
    'bool': get_bool,
    'int': get_int,
}
def tr(data, **kw): # pylint: disable=invalid-name
    """Translate between two dictionaries, apply a filter if needed"""
    for (source, dest) in kw.items():
        if source in data:
            value = data.pop(source)
            if isinstance(dest, tuple):
                # Translate the value using a predefined tr method
                value = TR_METHODS.get(dest[1], dest[1])(value, *dest[2:])
                dest = dest[0]
            if dest is not None and value not in ['', u'', None, 0]:
                data[dest] = value


def long_match(MAP, d, value, model=None, default='NOP', *cols, **filt):
    """Match in a model with case-insensitive multi-column matching."""
    queryset = filt.pop('queryset', None)
    if queryset is None and model is not None:
        queryset = model.objects
    match = filt.pop('_match', 'iexact')
    value = MAP.get(value, value)
    if value not in d and queryset is not None and cols:
        query = reduce(or_, [Q(**{'{}__{}'.format(col, match): value}) for col in cols])
        try:
            d[value] = queryset.filter(**filt).get(query)
        except queryset.model.DoesNotExist:
            if default != 'NOP':
                return default
            raise
    return d[value]


class StatusBar():
    """A generic command line status bar, use like so:

    for item in StatusBar("Label:", 300, looper_function):
        # Do something here.

    It will catch errors and return before the error prints.
    """
    io = sys.stderr # pylint: disable=invalid-name
    width = 40

    def __init__(self, label='', size=100, iterator=()):
        self.label = label
        self.size = size
        self.iter = iterator
        self.pos = 0
        self.down = self.width
        self.write("%s X%s]\r%s [" % (label, " " * self.width, label))

    def write(self, msg):
        """write the message to the status bar"""
        self.io.write(msg)
        self.io.flush()

    def count(self, item):
        """Return the count"""
        if self.count:
            return self.pos + 1
        return self.pos + len(item)

    def __iter__(self):
        try:
            for item in self.iter:
                self.pos = self.count(item)
                pon = int((self.pos / float(self.size)) * self.width)
                if self.width - pon < self.down:
                    self.write('-' * (self.down - self.width + pon))
                    self.down = self.width - pon
                yield item
        finally:
            self.write('X' * self.down + "\n")
