
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

