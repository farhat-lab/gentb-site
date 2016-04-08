
"""
Basic view mixins for predict views
"""

from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required

from .models import PredictDataset

class PredictMixin(object):
    """The baseline predict view"""
    slug_field = 'md5'
    model = PredictDataset

    @method_decorator(login_required)
    def dispatch(self, request, *args, **kwargs):
        """Only allow a logged in users to view"""
        return super(PredictMixin, self).dispatch(request, *args, **kwargs)

