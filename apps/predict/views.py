from __future__ import print_function

import os
import json

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext
from django.contrib.auth.decorators import login_required

from apps.predict.forms import UploadPredictionDataForm, SimpleConfirmationForm
from apps.predict.models import PredictDataset, PredictDatasetStatus,\
                                PredictDatasetNote, DatasetScriptRun
from apps.dropbox_helper.models import DropboxRetrievalLog
from apps.utils.view_util import get_common_dict
from apps.predict.message_helper import send_new_dataset_message_to_tb_admins
from django.views.generic import DetailView

from apps.utils.result_file_info import RESULT_OUTPUT_DIRECTORY_NAME, MATRIX_JSON_FILE_NAME

import logging
LOGGER = logging.getLogger(__name__)

#LoginRequiredMixin
class Heatmap(DetailView):
    queryset = PredictDataset.objects.filter(has_prediction=True)
    slug_field = 'md5'
    slug_url_kwarg = 'dataset_md5'
    template_name = 'predict/heatmap.html'


@login_required
def view_single_dataset(request, dataset_md5):

    d = get_common_dict(request, 'My Dataset Detail', view_my_datasets=True)

    # Not logged in, show login message
    #
    if not request.user.is_authenticated():
        return render_to_response('predict/my_dataset_detail.html',
                             d,
                             context_instance=RequestContext(request))

    try:
        dataset = PredictDataset.objects.get(md5=dataset_md5)
    except PredictDataset.DoesNotExist:
        raise Http404('PredictDataset not found')

    d['dataset'] = dataset
    d['dataset_notes'] = PredictDatasetNote.objects.filter(dataset=dataset).all()
    d['script_runs'] = DatasetScriptRun.objects.filter(dataset=dataset).all()
    d['tb_user'] = dataset.user

    # Pipeline command
    d['pipeline_command_found'] = dataset.get_pipeline_command()[0]

    try:
        local_directory_contents = os.listdir(dataset.file_directory)
        d['local_directory_contents'] = local_directory_contents
    except:
        d['local_directory_err'] = "Could not list directory contents"

    return render_to_response('predict/my_dataset_detail.html',
                             d,
                             context_instance=RequestContext(request))
@login_required
def view_my_datasets(request):
    d = get_common_dict(request, 'My Datasets', view_my_datasets=True)

    # Not logged in, show login message
    #
    if not request.user.is_authenticated():
        return render_to_response('predict/my_datasets.html',
                             d,
                             context_instance=RequestContext(request))


    d['datasets'] = PredictDataset.objects.filter(user=request.user.tbuser)

    return render_to_response('predict/my_datasets.html',
                             d,
                             context_instance=RequestContext(request))
