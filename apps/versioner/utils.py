#
# Copyright 2013, Martin Owens
#
# inkscape-web is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# inkscape-web is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with inkscape-web.  If not, see <http://www.gnu.org/licenses/>.
#

from django.template.context import Context

def to(t=list):
    """Create an object from a generator function, default is list"""
    def _outer(f):
        def _inner(*args, **kwargs):
            return t(f(*args, **kwargs))
        return _inner
    return _outer


def context_items(context):
    """Unpack a django context, equiv of dict.items()"""
    if not isinstance(context, Context):
        context = [context]
    for d in context:
        for (key, value) in d.items():
            yield (key, value)

class BaseMiddleware():
    """ 
    When used by a middleware class, provides a predictable get()
    function which will provide the first available variable from
    first the context_data, then the view, then the middleware.
    """
    def __init__(self, handler):
        self.get_response = handler

    def __call__(self, request):
        return self.get_response(request)

    @staticmethod
    def get(data, key, default=None, then=None):
        """Returns a data key from the context_data, the view, a get
        method on the view or a get method on the middleware in that order.
        
        Returns default (None) if all fail."""
        if key in data:
            return data[key]
        view = data.get('view', None)
        if hasattr(view, key):
            return getattr(view, key)
        if hasattr(view, 'get_'+key):
            return getattr(view, 'get_'+key)()
        if hasattr(then, 'get_'+key):
            return getattr(then, 'get_'+key)(data)
        return default

