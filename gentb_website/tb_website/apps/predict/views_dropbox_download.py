from __future__ import print_function
import os
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext

from apps.predict.forms import DropboxDownloadAttemptForm
from apps.predict.models import PredictDataset, PredictDatasetStatus,\
                                PredictDatasetNote, DatasetScriptRun,\
                                DropboxRetrievalLog
#from apps.shared_data.process_file_helper import get_process_file_results
#from apps.utils.view_util import get_common_dict
from apps.predict.message_helper import send_new_dataset_message_to_tb_admins



def record_file_retrieval_results(request):
    """
    Record the locations of the drobbox file retrieval attempt
    """
    if not request.POST:
        return HttpResponse('Method Not Allowed (use POST)', status=405)

    # --------------------------------------
    # Validate the results
    # --------------------------------------
    f = DropboxDownloadAttemptForm(request.POST)
    if not f.is_valid():
        data = json.dumps(dict(success=False,
                            message=f.errors))
        return HttpResponse(data, content_type='application/json')

    # --------------------------------------
    # Get the dataset
    # --------------------------------------
    try:
        dataset = PredictDataset.objects.get(md5=f.get_run_md5())
    except PredictDataset.DoesNotExist:
        error_msg = "A PredictDataset with was not found for md5: %s" % f.get_run_md5()
        data = json.dumps(dict(success=False,
                            message=error_msg)
                          )
        return HttpResponse(data, content_type='application/json')
