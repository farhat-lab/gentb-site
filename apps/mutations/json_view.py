
from django.views.generic import View
from django.core.serializers.json import DjangoJSONEncoder
from django.http import JsonResponse, HttpResponse

class DjangoJSONEncoder2(DjangoJSONEncoder):
    def default(self, obj):
        if isinstance(obj, timedelta):
            ARGS = ('days', 'seconds', 'microseconds')
            return {'__type__': 'datetime.timedelta',
                    'args': [getattr(obj, a) for a in ARGS]}
	return DjangoJSONEncoder.default(self, obj)

class JsonView(View):
    def dispatch(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        return self.render_to_response(context)

    def render_to_response(self, context, **kwargs):
        return JsonResponse(context, encoder=DjangoJSONEncoder2)

