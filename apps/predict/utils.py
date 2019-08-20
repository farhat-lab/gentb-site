# -*- coding: utf-8 -*-

from django.contrib.staticfiles.templatetags.staticfiles import static as _static
from django.utils.functional import lazy

def static(path, site_id=None):
    return _static(path)

def get_absolute_url_for_site(url, site):
    url_tmpl = "{scheme}//{domain}{url}"
    scheme = site.scheme and "{0}:".format(site.scheme) or ""
    return url_tmpl.format(scheme=scheme, domain=site.domain, url=url)

static_lazy = lazy(static, str)
