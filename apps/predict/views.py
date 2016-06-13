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

from .models import PredictDataset, DatasetScriptRun
from .mixins import PredictMixin, CallbackMixin
from .message_helper import send_dataset_run_message_to_tb_admins_and_user
from .forms import *


class Datasets(PredictMixin, ListView):
    pass


class DatasetView(PredictMixin, DetailView):
    pass


class Heatmap(PredictMixin, DetailView):
    queryset = PredictDataset.objects.filter(has_prediction=True)
    template_name = 'predict/heatmap.html'


class UploadChoices(PredictMixin, TemplateView):
    template_name = 'predict/predictdataset_upload.html'


class UploadManual(PredictMixin, CreateView):
    template_name = 'predict/predictdataset_manual.html'
    form_class = ManualInputForm

    def get_initial(self):
        return {
          'user': self.request.user,
          'status': PredictDataset.STATUS_FILE_RETRIEVAL_COMPLETE,
          'file_type': 'manual',
        }


class UploadView(PredictMixin, CreateView):
    @property
    def form_class(self):
        if self.kwargs['type'] == 'fastq':
            if self.kwargs['fastq'] == 'pair-end':
                return UploadFastQPairForm
            return UploadFastQSingleForm
        return UploadVcfForm

    def get_initial(self):
        return {
          'user': self.request.user,
          'status': PredictDataset.STATUS_CONFIRMED,
          'file_type': self.kwargs['type'],
          'fastq_type': self.kwargs.get('fastq', None),
        }


class Callback(CallbackMixin, FormView):
    """
    When the files have been processed, we callback and update statuses.
    """
    form_class = NotificationForm

    def failure(self, msg, status=400):
        return self.render_to_response({'success': False, 'message': msg},
                status=status)

    def form_valid(self, form):
        try:
            dataset_run = DatasetScriptRun.objects.get(md5=self.kwargs['slug'])
        except DatasetScriptRun.DoesNotExist:
            return self.failure("The ScriptRun was not found", 404)

        if dataset_run.result_received:
            return self.failure("This script has already run")

        dataset_run.result_received = True
        dataset_run.result_success = form.was_run_successful()
        dataset_run.result_data = form.get_result_data()
        dataset_run.save()

        send_dataset_run_message_to_tb_admins_and_user(dataset_run)
        #send_dataset_run_message_to_tb_admins(dataset_run)

        return self.render_to_response({'success': True, 'message': "OK"})

