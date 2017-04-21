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
Mixins specially for the maps app
"""

from datetime import timedelta

from django.core.serializers.json import DjangoJSONEncoder
from django.db.models.query import ValuesQuerySet
from django.views.decorators.cache import cache_page
from django.views.generic import View

from django.http import JsonResponse, HttpResponse
from django.conf import settings

class DjangoJSONEncoder2(DjangoJSONEncoder):
    """A json encoder to deal with the python objects we may want to encode"""
    def default(self, obj):
        if isinstance(obj, timedelta):
            ARGS = ('days', 'seconds', 'microseconds')
            return {'__type__': 'datetime.timedelta',
                    'args': [getattr(obj, a) for a in ARGS]}
        if isinstance(obj, ValuesQuerySet):
            return [item for item in obj]
        return DjangoJSONEncoder.default(self, obj)


class JsonView(View):
    """Quickly serve a python data structure as json"""
    cache_timeout = 5 * 60 * 60 * 24

    def get_cache_timeout(self):
        return self.cache_timeout

    def dispatch(self, *args, **kwargs):
        def _dispatch(request, *args, **kwargs):
            context = self.get_context_data(**kwargs)
            return self.render_to_response(context)
        return cache_page(self.get_cache_timeout())(_dispatch)(*args, **kwargs)

    def render_to_response(self, context, **kwargs):
        return JsonResponse(context, encoder=DjangoJSONEncoder2)


class DataSlicerMixin(object):
    """
    Provide a way to slice up a given model based on inputs.
    """
    filters = {}
    required = []
    order = []
    values = []

    def get_model(self):
        """Return the basic model to slice"""
        try:
            return self.model
        except:
            raise NotImplementedError("You must provide a model to slice.")

    def get_filters(self, without=None):
        """Gets the filter applied to the queryset for this slice."""
        for (key, filtr) in self.filters.items():
            if without and filtr.startswith(without):
                continue
            if key not in self.request.GET and key not in self.required:
                continue
            if filtr.endswith('__in'):
                yield (filtr, self.request.GET.getlist(key, []))
            else:
                yield (filtr, self.request.GET.get(key, ''))

    def get_queryset(self, without=None):
        """Applies any filters from the request query to the given model"""
        qs = self.get_model().objects.all()
        if self.order:
            qs = qs.order_by(*self.order)
        return qs.filter(**dict(self.get_filters(without)))

    def get_data(self, without=None):
        qs = self.get_queryset(without=without)
        if self.values:
            vals = [v for v in self.values if v != without]
            qs = qs.values(*vals)
        return qs

    def get_list(self, qs, column):
        """Returns a flat list for this column"""
        return qs.values_list(column, flat=True).distinct().order_by(column)

