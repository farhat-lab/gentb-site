# -*- coding: utf-8 -*-

from django.contrib.staticfiles.templatetags.staticfiles import static as _static
from django.utils.functional import lazy

try:
    # For django >= 2.0
    from django.urls import reverse as _reverse
except ImportError:
    from django.core.urlresolvers import reverse as _reverse

def static(path, site_id=None):
    return _static(path)

def get_absolute_url_for_site(url, site):
    url_tmpl = "{scheme}//{domain}{url}"
    scheme = site.scheme and "{0}:".format(site.scheme) or ""
    return url_tmpl.format(scheme=scheme, domain=site.domain, url=url)

static_lazy = lazy(static, str)
