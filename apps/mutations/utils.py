#
# Copyright (C) 2016   Dr. Maha Farhat
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
import json
import sys
import time

from collections import defaultdict, OrderedDict
from datetime import date

MONTHS = ['', 'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
DATE_FORMATS = [
   r'^(?P<year>\d{4})-(?P<mm>\d{2})-(?P<dd>\d{2})$',
   r'^(?P<year>\d{2,4})$',
   r'^(?P<mm>\d{1,2})/(?P<year>\d{2,4})$',
   r'^(?P<mon>[A-Za-z]{3})/(?P<year>\d{2,4})$',
   r'^(?P<mm>\d{1,2})/(?P<dd>\d{1,2})/(?P<year>\d{2,4})$',
   r'^(?P<dd>\d{2})(?P<mon>[A-Za-z]{3})(?P<year>\d{2,4})$',
]

def re_match(re_list, string, raw=False):
    """
    Match the first regular expression in a list of regular expressions
    and return the groupdict of the matching expression.

    Will return raw match object and index in list if raw is True.

    Will replace your list items with re objects in place! for caching.
    """
    string = unicode(string).strip()
    for (x, r) in enumerate(re_list):
        # Compile only once per run
        if isinstance(r, str):
            r = re.compile(r)
            re_list[x] = r
        match = r.search(string)
        if match is not None:
            return x, match if raw else match.groupdict()
    raise ValueError("Couldn't match string: %s" % str(string))

def re_match_dict(re_dict, string):
    """
    Like re_match(), but instead of a list, it takes a dictionary. The first
    matching regular expression then has it's value populated with the
    match's group dictionary as a template. So for example:
    
    re_match_dict({r'P(\d+)': 'A%s'}, 'P45') == 'A45'

    It will also take a list of tuples for ordered matching.
    """
    string = unicode(string).strip()
    re_dict = re_dict.items() if isinstance(re_dict, dict) else re_dict
    for (r, value) in re_dict:
        # Compile only once per run
        if isinstance(r, str):
            r = re.compile(r)
        match = r.search(string)
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

mutation_re = [
      r'^SNP_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_(?P<coding>[ACTG]\d+[ACTG])_(((?P<amino>[A-Z\*]\d+[A-Z\*])|(?P<noncode>promoter|inter))_(?P<gene>[a-zA-Z\d\-_]+)|(?P<rgene>rr[sl]))\'?$',
      r'^(?P<mode>(INS|DEL))_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_(i|d|\.|i\.)?(?P<codes>[\d\-]+[ATGC]*)_((?P<noncode>promoter|inter|\d+)_)?(?P<gene>[a-z][a-zA-Z\d\-_]+?)(_(?P<amino>[A-Z\*]\d+[A-Z\*]))?\'?$',
      # These are older SNP names and should probably be converted
      r'SNP_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_(?P<coding>[ACTG]\d+[ACTG])_(?P<gene>Rv\w+)',
      r'SNP_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_(?P<coding>[ACTG]\d*[ACTG])_PE_(?P<gene>[a-zA-Z\d\-]+)_(?P<amino>[A-Z\*]\d+[A-Z\*])',
      r'SNP_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_(?P<coding>[ACTG]\d*[ACTG])_(?P<gene>[a-zA-Z\d\-]+)_(?P<amino>[A-Z\*]\d+[A-Z\*])',
      r'SNP_(?P<syn>[A-Z]{1,2})_(?P<ntpos>\d+)_(?P<coding>[ACTG]\d*[ACTG])_(?P<gene>[a-zA-Z\d\-_]+)',
]

def match_snp_name(name, **kw):
    """
    Tries to regex match the snp name and returns a dictionary of information decoded from the SNP name.
    """
    return re_match(mutation_re, name, **kw)

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
    index, match = match_snp_name(mutation, raw=True)
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

    return ''.join(ret), mutation_re[index].pattern, snp


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

    if snp['syn'] == 'P':
	if snp.get('noncode', '') != 'promoter':
	    _ = ValueError("Promoter doesn't specify 'promoter' part in %s" % name)
	return (index, 'promoter ' + gene, name)
    elif snp['syn'] == 'I':
	if snp.get('noncode', '') != 'inter':
	    _ = ValueError("Integenic doesn't specify 'inter' part")
	return (index, 'intergenic ' + gene, name)
    elif snp['syn'] in ['CN', 'CD', 'CF', 'CI', 'CS', 'CZ', 'N', 'ND', 'NI', 'NF']:
	return (index, gene, name)

    raise ValueError("Must be promoter, intergenic or CN, CD, CF, CI, CS, CZ or N, ND, NI, NF")


class defaultlist(defaultdict):
    """Like a defaultdict, but generates a list of items with the same key"""
    def __init__(self, generator):
        super(defaultlist, self).__init__(list)
	for key, value in generator:
	    self[key].append(value)

    def flatten(self, method, **kw):
        """Trun each list value into a single value based on method"""
        for key, value in self.items():
            self[key] = method(value, **kw)

    def re_key(self, new_key, pop=True):
        """Assuming each value is a dictionary, re-key the dictionary"""
        for key in list(self):
            value = self.pop(key)
            self[value.pop(new_key)].append(value)


def json_generator(f):
    """Decorate a method that processes on row in a json filename list"""
    def _inner(*args, **kw):
        args = list(args)
        index = 0
        if not isinstance(args[0], str) and len(args) > 1:
            index = 1
        if not os.path.isfile(args[index]):
            raise IOError("Json file '%s' Not Found" % args[index])
        with open(args[index], 'r') as fhl:
            rows = json.loads(fhl.read())
            if 'status' in kw:
                rows = StatusBar(kw.pop('status'), len(rows), rows, True)
            for row in rows:
                args[index] = row
                value = f(*args, **kw)
                if value is not None:
                    yield value
    return _inner

def to(method):
    """Turn generators into objects, method can be a type, object or function"""
    def _outer(f):
        def _inner(*args, **kw):
            return method(f(*args, **kw))
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
    return str(string).upper() not in ['N', 'NO', 'FALSE', 'F', '0', 'NONE']

def get_int(datum):
    if isinstance(datum, basestring):
        try:
            return int(datum)
        except ValueError:
            if len(datum) == 1:
                return ord(datum)
    return datum

tr_methods = {
  'date': get_date,
  'bool': get_bool,
  'int': get_int,
}
def tr(data, **kw):
    """Translate between two dictionaries, apply a filter if needed"""
    for (source, dest) in kw.items():
        if source in data:
            value = data.pop(source)
            if isinstance(dest, tuple):
                # Translate the value using a predefined tr method
                value = tr_methods.get(dest[1], dest[1])(value, *dest[2:])
                dest = dest[0]
            if dest is not None and value not in ['', u'', None, 0]:
                data[dest] = value


class StatusBar(object):
    """A generic command line status bar, use like so:

    for item in StatusBar("Label:", 300, looper_function):
        # Do something here.

    It will catch errors and return before the error prints.
    """
    io = sys.stderr
    width = 40

    def __init__(self, label='', size=100, iterator=[], count=False):
        self.label = label
        self.size = size
        self.iter = iterator
        self.pos = 0 
        self.down = self.width
        self.write("%s X%s]\r%s [" % (label, " " * self.width, label))

    def write(self, msg):
        self.io.write(msg)
        self.io.flush()

    def count(self, item):
        if self.count:
            return self.pos + 1
        return self.pos + len(item)

    def __iter__(self):
        try:
            for item in self.iter:
                self.pos = self.count(item)
                p = int((self.pos / float(self.size)) * self.width)
                if self.width - p < self.down:
                    self.write('-' * (self.down - self.width + p)) 
                    self.down = self.width - p 
                yield item
        finally:
            self.write('X' * self.down + "\n")

