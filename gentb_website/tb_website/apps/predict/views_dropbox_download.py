from __future__ import print_function
import os
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext

from apps.predict.forms import DropboxDownloadAttemptForm
from apps.predict.models import PredictDataset, PredictDatasetStatus,\
                                PredictDatasetNote, DatasetScriptRun
from apps.dropbox_helper.models import DropboxRetrievalLog
from apps.predict.message_helper import send_new_dataset_message_to_tb_admins



def record_file_retrieval_results(request):
    """
    Record the locations of the dropbox file retrieval attempt
    """
    if not request.POST:
        return HttpResponse('Method Not Allowed (use POST)',
            status=405,
            content_type='application/json')

    # --------------------------------------
    # Validate the results
    # --------------------------------------
    f = DropboxDownloadAttemptForm(request.POST)
    if not f.is_valid():
        data = json.dumps(dict(success=False,
                            message=f.errors))
        return HttpResponse(data, content_type='application/json')

    # --------------------------------------
    # Retrieve the DropboxRetrievalLog
    # --------------------------------------
    dbox_log = DropboxRetrievalLog.objects.filter(md5=f.get_run_md5()).first()
    if dbox_log is None:
        error_msg = "A DropboxRetrievalLog was not found. (md5 of log: %s" % dataset.md5
        data = json.dumps(dict(success=False,
                            message=error_msg)
                          )
        return HttpResponse(data, content_type='application/json')

    # --------------------------------------
    # Log the results
    # --------------------------------------
    dbox_log.set_retrieval_end_time()
    if not f.was_run_successful():     # Failure!
        success_flag = False
        msg = "There was an error retrieving the dropbox files."
        dbox_log.dataset.set_status_file_retrieval_error()
        dbox_log.files_retrieved = False
        dbox_log.retrieval_error = f.get_result_data()
    else:                               # Success!
        success_flag = True
        msg = "Dropbox files successfully retrieved."
        dbox_log.dataset.set_status_file_retrieval_complete()
        dbox_log.files_retrieved = True

    dbox_log.save()

    data = json.dumps(dict(success=success_flag,
                    message=msg))
    return HttpResponse(data, content_type='application/json')
