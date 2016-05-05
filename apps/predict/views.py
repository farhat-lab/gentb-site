"""
Views for predict app
"""
from __future__ import print_function

import logging

LOGGER = logging.getLogger(__name__)

from django.views.generic import (
    DetailView, ListView, CreateView, UpdateView, FormView
)
from django.core.urlresolvers import reverse

from .models import PredictDataset, DatasetScriptRun
from .forms import UploadForm, UploadConfirmForm, NotificationForm
from .mixins import PredictMixin, CallbackMixin
from .message_helper import send_dataset_run_message_to_tb_admins_and_user

class Datasets(PredictMixin, ListView):
    pass

class DatasetView(PredictMixin, DetailView):
    pass

class Heatmap(PredictMixin, DetailView):
    queryset = PredictDataset.objects.filter(has_prediction=True)
    template_name = 'predict/heatmap.html'


class UploadView(PredictMixin, CreateView):
    form_class = UploadForm

    def get_initial(self):
        return {'status': PredictDataset.STATUS_NOT_READY, 'user': self.request.user}

    def get_success_url(self):
        return reverse('predict:view_predict_upload_step2_confirm', kwargs=dict(slug=self.object.md5))


class UploadConfirm(PredictMixin, UpdateView):
    """
    After dropbox has downloaded some metadata, we confirm the files.
    """
    template_name = 'predict/predictdataset_confirm.html'
    form_class = UploadConfirmForm


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

