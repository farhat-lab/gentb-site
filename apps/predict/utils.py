# -*- coding: utf-8 -*-

import csv
from io import StringIO

from django.contrib.staticfiles.templatetags.staticfiles import static as _static
from django.utils.functional import lazy

def static(path, site_id=None):
    return _static(path)
static_lazy = lazy(static, str)

def get_absolute_url_for_site(url, site):
    url_tmpl = "{scheme}//{domain}{url}"
    scheme = site.scheme and "{0}:".format(site.scheme) or ""
    return url_tmpl.format(scheme=scheme, domain=site.domain, url=url)

def lineage_spoligo(data):
    return [('spoligo', data[1] + " / " + data[-1])]

def lineage_fast_caller(data):
    reader = csv.reader(StringIO(data), delimiter="\t")
    for name, value in zip(next(reader), next(reader)):
        if 'Isolate' in name:
            continue
        yield (name, value)

def lineage_other_caller(data):
    data = [lin.replace('lineage', '') for lin in data.split(',')]
    for x, lin in enumerate(data):
        if x and not lin.startswith(data[x-1]):
            return [('lineage', data)]
    return [('lineage', data[-1])]

def filter_none(vals):
    """Remove none values"""
    ret = []
    for x in vals:
        if x in (u'None', u'Null'):
            x = None
        ret.append(x)
    return ret
