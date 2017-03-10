"""
Views for predict app
"""
from __future__ import print_function

import logging

LOGGER = logging.getLogger(__name__)

from django.views.generic import (
    DetailView, ListView, CreateView, UpdateView, FormView,
    TemplateView,
)
from django.core.urlresolvers import reverse
from django.http.response import JsonResponse
from django.http import Http404

from .models import PredictDataset, PredictDatasetNote
from .mixins import PredictMixin
from .forms import *


class Datasets(PredictMixin, ListView):
    title = "My Datasets"

    @classmethod
    def get_absolute_url(cls):
        return reverse('predict:view_my_datasets')


class DatasetView(PredictMixin, DetailView):
    parent = Datasets

    def get_object(self):
        obj = super(DatasetView, self).get_object()
        for strain in obj.strains.all():
            strain.update_status()
        return obj

class Heatmap(PredictMixin, DetailView):
    queryset = PredictDataset.objects.all()
    template_name = 'predict/heatmap.html'


class UploadChoices(PredictMixin, TemplateView):
    template_name = 'predict/predictdataset_upload.html'
    forms = UploadForm.all_forms()
    title = "Create Prediction"
    parent = Datasets

    @classmethod
    def get_absolute_url(cls):
        return reverse('predict:upload')


class UploadView(PredictMixin, CreateView):
    model = PredictDataset
    parent = UploadChoices

    def get_title(self):
        return self.form_class.title

    def get_template_names(self):
        default = super(UploadView, self).get_template_names()
        default = ['predict/predictdataset_%s.html' % self.kwargs['type']] + default
        return default

    @property
    def form_class(self):
        for form in UploadForm.all_forms():
            if form.my_file_type == self.kwargs['type']:
                return form
        raise Http404("No input type: %s" % self.kwargs['type'])

    def get_initial(self):
        return {
          'user': self.request.user,
          'status': self.form_class().my_status,
          'file_type': self.kwargs['type'],
        }


class AddNote(PredictMixin, CreateView):
    model = PredictDatasetNote
    fields = ('note',)

    def form_invalid(self, form):
        return JsonResponse({'status': 'INVALID'})

    def form_valid(self, form):
        obj = form.save(commit=False)
        obj.title = unicode(self.request.user)
        obj.dataset = self.get_queryset().get(md5=self.kwargs['slug'])
        obj.save()
        return JsonResponse({
          'status': 'OK',
          'title': obj.title,
          'note': obj.note,
        })

