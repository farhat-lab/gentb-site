"""
Shows all the available urls for a django website, useful for debugging.
"""

import types

from django.core.management.base import BaseCommand

import tb_website.urls

class Url():
    is_module = False
    is_view = False

    def __init__(self, parent, entry, module):
        self.depth = 0 
        self.parent = parent
        self.entry = entry
        self.pattern = entry.regex.pattern
        self.module = module

        if parent is not None:
            self.depth = parent.depth + 1 

    def __str__(self):
        if self.name:
            return "%s [%s]" % (self.full_pattern, self.name)
        return self.full_pattern

    @property
    def name(self):
        name = getattr(self.entry, 'name', None)
        namespace = self.namespace
        if name and namespace:
            return "%s:%s" % (namespace, name)
        if name:
            return name
        return None

    @property
    def slug(self):
        name = self.name
        if name is None:
            return slugify(self.pattern)
        if self.kwargs:
            name += '?' + '+'.join(self.kwargs)
        return name

    @property
    def kwargs(self):
        """Gathers all kwargs from this and every parent regex"""
        kw = self.parent.kwargs if self.parent else {}
        kw.update(self.entry.regex.groupindex)
        return kw

    @property
    def full_pattern(self):
        pattern = self.pattern.lstrip('^').rstrip('$')
        if self.parent:
            if pattern and self.parent.pattern.endswith('$'):
                logger.warning("Possible broken url, parent url ends "
                        "string matching: " + str(self.parent))
            pattern = self.parent.full_pattern + pattern
        if pattern == 'None/':
            return '/'
        return pattern

    @property
    def namespace(self):
        if hasattr(self.entry, 'namespace') and self.entry.namespace:
            return self.entry.namespace
        if self.parent is not None:
            return self.parent.namespace
        return None

    def test_url(self, *args, **kw):
        """Atttempt to generate a test url based on these kwargs"""
        if self.name:
            return reverse(self.name, args=args, kwargs=kw)
        if not self.kwargs:
            # Construct the pattern without any
            return '/' + self.full_pattern.lstrip('^').rstrip('$')
        return None

class UrlModule(Url):
    """A url include"""
    is_module = True

    def __str__(self):
        tag = super(UrlModule, self).__str__()
        return "%s > %s >>" % (self.full_pattern, self.name)

    @property
    def name(self):
        return self.urls_name(self.module)

    def urls_name(self, uc):
        if isinstance(uc, list) and uc:
            return self.urls_name(uc[0])
        if hasattr(uc, '__name__'):
            return uc.__name__
        return None

class UrlFunction(Url):
    """A url using a simple function"""
    def __str__(self):
        tag = super(UrlFunction, self).__str__()
        return "%s > %s()" % (tag, self.module.__name__)

class WebsiteUrls():
    """A class that can loop through urls in a tree structure"""
    def __iter__(self):
        """
        Yields every url with a Url class, see Url() for details.
        """
        dupes = {}
        for item in self.url_iter(tb_website.urls.urlpatterns):
            key = (item.name, item.full_pattern)
            if not item.is_module:
                if key is not None and key in dupes:
                    logger.error(
                       "URL Name is already used '%s' -> '%s'" % key\
                       + "\n  a) " + str(dupes[key])\
                       + "\n  b) " + str(item) + '\n')
                dupes[key] = item
            yield item

    def url_iter(self, urllist, parent=None):
        """
        Returns a specific arm of the urls tree (see __iter__)
        """
        for entry in urllist:
            if hasattr(entry, 'url_patterns'):
                if hasattr(entry, '_urlconf_module'):
                    this_parent = UrlModule(parent, entry, entry._urlconf_module)
                    if this_parent.name:
                        yield this_parent

                    # replace with yield from (python 3.3) when possible.
                    for item in self.url_iter(entry.url_patterns, this_parent):
                        yield item

                continue

            if hasattr(entry, '_callback'):
                callback = entry._callback
                if isinstance(callback, types.FunctionType):
                    yield UrlFunction(parent, entry, callback)
                elif hasattr(callback, 'model') and callback.model is not None:
                    yield UrlView(parent, entry, callback)
                else:
                    yield Url(parent, entry, callback)


class Command(BaseCommand):
    args = '<start_url>'
    help = 'Shows all urls begining with the start_url'

    def handle(self, *args, **options):
        self.start_url = None
        if len(args) > 0:
            self.start_url = args[0]

        for url in WebsiteUrls():
            self.stdout.write(" " * url.depth + str(url))

