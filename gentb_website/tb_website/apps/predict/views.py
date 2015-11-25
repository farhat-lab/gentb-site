from __future__ import print_function
import os
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
from django.contrib.auth.decorators import login_required
from apps.predict.pipeline_hardcoded_script_runner import PipelineScriptRunner

import logging
logger = logging.getLogger(__name__)

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
    (pipeline_command_found, pipeline_command_or_err) = PipelineScriptRunner.get_pipeline_command(dataset)
    d['pipeline_command_found'] = pipeline_command_found

    if pipeline_command_found:
        d['pipeline_command'] = pipeline_command_or_err
    else:
        d['pipeline_command_err'] = pipeline_command_or_err

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
    """
    # has the script already been run successfully?
    if shared_file_info.prediction_results:
        # Yes, use those results
        #
        msg_or_data = shared_file_info.get_prediction_results()
        if msg_or_data:
            success = True
        else:
            success = False
    else:
        # No, then run the prediction R script
        #
        (success, msg_or_data) = get_process_file_results(shared_file_info.file_obj.file.name)
        if success is True:
            shared_file_info.set_prediction_results(msg_or_data)
            shared_file_info.save()

    #print ('msg_or_data', type(msg_or_data))
    d['UPLOAD_SUCCESS'] = True

    d['FILE_PROCESS_SUCCESS'] = success
    d['FILE_PROCESS_ERR_OR_DATA'] = msg_or_data

    d['shared_file_info'] = shared_file_info


    return render_to_response('predict/predict_upload_success.html',
                             d,
                             context_instance=RequestContext(request))
    """
