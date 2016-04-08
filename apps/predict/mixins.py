"""
Basic view mixins for predict views
"""

from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse

from .remote_forms import RemoteForm
from .models import PredictDataset

class PredictMixin(object):
    """The baseline predict view"""
    slug_field = 'md5'
    model = PredictDataset

    @method_decorator(login_required)
    def dispatch(self, request, *args, **kwargs):
        """Only allow a logged in users to view"""
        return super(PredictMixin, self).dispatch(request, *args, **kwargs)


class CallbackMixin(object):
    """A view called by the server itself"""
    slug_field = 'md5'
    model = PredictDataset

    @method_decorator(csrf_exempt)
    def dispatch(self, request, *args, **kwargs):
        """Test IP address instead of csrf token"""
        # XXX - Test ip address of requestor here
        return super(CallbackMixin, self).dispatch(request, *args, **kwargs)

    def render_to_response(self, context, **kw):
        """Return a json object instead of a template"""
        context.pop('view', None)
        context['form'] = RemoteForm(context['form']).as_dict()
        return JsonResponse(context, **kw)

