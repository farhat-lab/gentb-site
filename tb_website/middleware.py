#
# Copyright 2016, Martin Owens <doctormo@gmail.com>
#           2017, Maha Farhat for GenTB
#
# This file is part of the software inkscape-web, consisting of custom 
# code for the Inkscape project's django-based website.
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
# Copied from inkscape-web 2017-01-20 and modified for GenTB
#
"""
Extra functionality for getting titles and breadcrumbs from views.
"""

from django.urls import reverse
from django.utils.translation import ugettext_lazy as _

from django.db.models import Manager, QuerySet
from django.contrib.contenttypes.models import ContentType

def IterObject(typ=list): # pylint: disable=invalid-name
    """Create an object from a generator function, default is list"""
    def __outer__(func):
        def __inner__(*args, **kwargs):
            return typ(func(*args, **kwargs))
        return __inner__
    return __outer__

class AutoBreadcrumbMiddleware():
    """
    This middleware controls and inserts some breadcrumbs
    into most pages. It attempts to navigate object hierachy
    to find the parent objects and their natural links.

    This can get kind of complicated so most things can be overloaded
    by the view that's being called.

    As well as filling in a 'breadcrumbs' variable, it will also fill in
    a 'title' variable for general use. Each view can contain a get_title
    method or title propery to easily insert a title.
    """
    keys = ('breadcrumbs', 'title')

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)
        return response

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

    def process_template_response(self, request, response):
        if not hasattr(response, 'context_data'):
            return response
        data = response.context_data
        for key in self.keys:
            response.context_data[key] = self.get(data, key, then=self)
        if len(response.context_data['breadcrumbs']) == 1:
            del response.context_data['breadcrumbs']
        response.context_data['admin_link'] = self.get_admin_link(request.user, data)
        return response

    @staticmethod
    def get_admin_link(user, data):
        """Generates an admin link to edit this object"""
        obj = data.get('object', None)
        if obj is not None and user is not None:
            ct = ContentType.objects.get_for_model(type(obj))
            args = (obj.pk,) if obj else ()
            if user.has_perm(f'{ct.app_label}.change_{ct.model}'):
                return {
                    'name': f'Edit "{obj}"',
                    'url': reverse(f'admin:{ct.app_label}_{ct.model}_change', args=args),
                }

    @staticmethod
    def get_title(data):
        """If no title specified in context, use last breadcrumb"""
        if data.get('breadcrumbs', False):
            return list(data['breadcrumbs'])[-1][-1]
        return None

    @IterObject(list)
    def get_breadcrumbs(self, data):
        """Return breadcrumbs only called if no breadcrumbs in context"""
        title = self.get(data, 'title')
        parent = self.get(data, 'parent')
        obj = self.get(data, 'object')
        page = self.get(data, 'current_page')
        if not obj and page:
            # django-cms pages already have Home
            obj = page
        else:
            yield (reverse('home'), _('Home'))

        root = self.get(data, 'breadcrumb_root')
        if root:
            if isinstance(root, list):
                for item in root:
                    yield self.object_link(item)
            else:
                yield self.object_link(root)

        lst = self.get(data, 'object_list')
        if isinstance(lst, (Manager, QuerySet)):
            if obj is None:
                pass # obj = lst
            elif parent is None:
                parent = lst

        for obj in self.get_ancestors(obj, parent):
            link = self.object_link(obj)
            if link is not None:
                yield link

        if title is not None:
            yield (None, title)

    def get_ancestors(self, obj, parent=None):
        if hasattr(obj, 'breadcrumb_parent'):
            parent = obj.breadcrumb_parent()
        else:
            parent = getattr(obj, 'parent', parent)

        if parent is not None:
            for ans in self.get_ancestors(parent):
                yield ans
        yield obj

    @staticmethod
    def object_link(obj):
        """Get name from object model"""
        url = None
        if obj is None or (isinstance(obj, tuple) and len(obj) == 2):
            return obj
        if hasattr(obj, 'breadcrumb_name'):
            name = obj.breadcrumb_name()
        elif hasattr(obj, 'name'):
            name = obj.name
        elif hasattr(obj, 'title'):
            name = obj.title
        else:
            try:
                name = str(obj, errors='ignore')
            except UnicodeEncodeError:
                name = "Name Error"
        if hasattr(obj, 'get_absolute_url'):
            url = obj.get_absolute_url()
        if name is not None and name.startswith('['):
            return None
        return (url, name)
