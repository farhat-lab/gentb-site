from __future__ import print_function
import os
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext
from django.contrib.auth.decorators import login_required

from apps.predict.forms import UploadPredictionDataForm, SimpleConfirmationForm
from apps.predict.models import PredictDataset, PredictDatasetStatus,\
                                PredictDatasetNote, DatasetScriptRun,\
                                DropboxRetrievalLog
from apps.shared_data.process_file_helper import get_process_file_results
from apps.utils.view_util import get_common_dict
from apps.predict.message_helper import send_new_dataset_message_to_tb_admins
from django.contrib.auth.decorators import login_required
from subprocess import Popen


def view_predict_page(request):

    d = get_common_dict(request, 'Predict', predict_page=True)

    # Not logged in, show login message
    #
    if not request.user.is_authenticated():
        return render_to_response('predict/predict_upload.html',
                             d,
                             context_instance=RequestContext(request))

    if request.POST:
        f = UploadPredictionDataForm(request.POST)
        if f.is_valid():
            new_dataset = f.get_dataset(request.user.tbuser)

            #kick_off_retrieval(new_dataset)
            #get_dropbox_metadata(new_dataset)
            #send_new_dataset_message_to_tb_admins(new_dataset)

            success_url = reverse('view_predict_upload_step2_confirm',
                                  kwargs=dict(dataset_md5=new_dataset.md5)
                                )
            return HttpResponseRedirect(success_url)
        else:
            d['ERROR_FOUND']  = True
    else:
        f = UploadPredictionDataForm(label_suffix='')

    d['predict_form'] = f

    return render_to_response('predict/predict_upload.html',
                             d,
                             context_instance=RequestContext(request))


@login_required
def view_predict_upload_step2_confirm(request, dataset_md5):
    """
    At this point:
        - We have a PredictDataset object and
        - Have validated the Dropbox link--it has either .fastq or .vcf files

    Now, Ask the user to confirm that these are the files
    """
    d = get_common_dict(request, 'Predict Upload - Confirm', predict_page=True)

    # Retrieve the PredictDataset
    #
    try:
        predict_dataset = PredictDataset.objects.get(md5=dataset_md5)
    except PredictDataset.DoesNotExist:
        raise Http404('PredictDataset not found')

    d['dataset'] = predict_dataset
    d['tb_user'] = predict_dataset.user


    # Retrieve the DropboxRetrievalLog
    #
    dbox_log = DropboxRetrievalLog.objects.filter(dataset=predict_dataset).first()
    if dbox_log is None:
        raise Http404('DropboxRetrievalLog for PredictDataset not found')

    # Handle this better!
    if dbox_log.file_metadata_err_msg:
        return HttpResponse('DropboxRetrievalLog shows fail: {0}<br />{1}'.format(\
                        dbox_log.file_metadata_err_msg, 'view_predict_upload_step2_confirm')
                        )

    # Show or validate the confirmation form
    #
    if request.POST:
        f = SimpleConfirmationForm(request.POST)
        if f.is_valid():
            if f.do_not_use_files():
                """
                Delete the info and return to predict Upload
                """
                predict_dataset.delete() # cascades to DropboxRetrievalLog
                next_url = reverse('view_predict_upload_delete', args=())
            else:
                # >> kick off the download process here...
                next_url = reverse('view_predict_upload_success',
                                  kwargs=dict(dataset_md5=predict_dataset.md5)
                                )
            return HttpResponseRedirect(next_url)
        else:
            d['ERROR_FOUND']  = True
    else:
        f = SimpleConfirmationForm()

    d['dbox_log'] = dbox_log
    d['confirm_form'] = f

    return render_to_response('predict/predict_upload_step2.html',
                        d,
                        context_instance=RequestContext(request))


@login_required
def view_predict_upload_delete(request):
    return HttpResponse('view_predict_upload_delete')

@login_required
def view_predict_upload_success(request, dataset_md5):

    d = get_common_dict(request, 'Predict Upload Success', predict_page=True)

    try:
        dataset = PredictDataset.objects.get(md5=dataset_md5)
    except PredictDataset.DoesNotExist:
        raise Http404('PredictDataset not found')

    d['dbox_log'] = DropboxRetrievalLog.objects.filter(dataset=dataset).first()
    print ('get it?', d['dbox_log'])

    d['dataset'] = dataset
    d['tb_user'] = dataset.user

    return render_to_response('predict/predict_upload_success.html',
                             d,
                             context_instance=RequestContext(request))


def view_test_upload_success(request):


    #print ('msg_or_data', type(msg_or_data))
    d = get_common_dict(request, 'Test Predict Upload Success', predict_page=True)

    d['UPLOAD_SUCCESS'] = True

    d['FILE_PROCESS_SUCCESS'] = True
    d['FILE_PROCESS_ERR_OR_DATA'] = { 'test_data' : True }#msg_or_data

    #d['shared_file_info'] = shared_file_info


    return render_to_response('predict/predict_upload_success.html',
                             d,
                             context_instance=RequestContext(request))
