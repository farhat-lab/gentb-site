"""
Views for predict app
"""
from __future__ import print_function

import logging

LOGGER = logging.getLogger(__name__)

from apps.predict.models import PredictDataset
from django.views.generic import DetailView, ListView, CreateView

from .forms import UploadPredictionDataForm
from .mixins import PredictMixin

class Datasets(PredictMixin, ListView):
    pass

class DatasetView(PredictMixin, DetailView):
    pass

class Heatmap(PredictMixin, DetailView):
    queryset = PredictDataset.objects.filter(has_prediction=True)
    template_name = 'predict/heatmap.html'


class UploadView(PredictMixin, CreateView):
    form_class = UploadPredictionDataForm

    def get_initial(self):
        return {'status': PredictDataset.NOT_READY, 'user': self.request.user}

