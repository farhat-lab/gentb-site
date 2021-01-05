"""
Basic view mixins for predict views
"""

from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required

from .models import PredictDataset

class PredictMixin():
    """The baseline predict view"""
    slug_field = 'md5'

    @method_decorator(login_required)
    def dispatch(self, request, *args, **kwargs):
        """Only allow a logged in users to view"""
        return super(PredictMixin, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        """Limit queryset to the user's own predictions only"""
        qset = PredictDataset.objects.all()
        if 'slug' not in self.kwargs and not self.request.user.is_staff:
            # Limit to my own predictions unless I have the md5
            qset = qset.filter(user_id=self.request.user.pk)
        return qset.prefetch_related('strains', 'strains__piperun', 'strains__piperun__programs')

