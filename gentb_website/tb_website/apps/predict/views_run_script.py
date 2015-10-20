from __future__ import print_function
import json

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext
from django.views.decorators.csrf import csrf_exempt

from apps.predict.models import PredictDataset, PredictDatasetStatus, PredictDatasetNote, DatasetScriptRun
from apps.utils.view_util import get_common_dict, IS_ACTIVE_STAFF
from apps.predict.script_run_helper import run_script_on_dataset
from apps.predict.forms import DatasetRunNotificationForm

from apps.predict.message_helper import send_dataset_run_message_to_tb_admins

def view_run_dataset_script(request, dataset_md5):

    d = get_common_dict(request, 'Run Dataset Script', predict_page=True)

    if not request.user.is_authenticated():
        raise Http404('throw a 404 if not logged in--should not see link')

    if d.get(IS_ACTIVE_STAFF, False) is False:
        return render_to_response('predict/view_run_dataset_script.html',
                             d,
                             context_instance=RequestContext(request))

    # Get the dataset
    #
    try:
        dataset = PredictDataset.objects.get(md5=dataset_md5)
    except PredictDataset.DoesNotExist:
        raise Http404('PredictDataset not found')

    d['dataset'] = dataset

    # create a run object
    #
    (success, err_or_script_run) = run_script_on_dataset(dataset)

    return HttpResponseRedirect(reverse('view_single_dataset', kwargs=dict(dataset_md5=dataset.md5)))

    #return render_to_response('predict/view_run_dataset_script.html',
    #                         d,
    #                         context_instance=RequestContext(request))


@csrf_exempt
def view_dataset_run_notification(request):
    """
    Should have middleware to limit IPs with access to method
    """

    if not request.POST:
        return HttpResponse('Method Not Allowed (use POST)', status=405)

    f = DatasetRunNotificationForm(request.POST)
    if not f.is_valid():
        data = json.dumps(dict(success=False,
                            message=f.errors))
        return HttpResponse(data, content_type='application/json')

    print (f.cleaned_data)
    # --------------------------------------
    # Get the dataset run
    # --------------------------------------
    try:
        dataset_run = DatasetScriptRun.objects.get(md5=f.get_run_md5())
    except DatasetScriptRun.DoesNotExist:
        error_msg = "A DatasetScriptRun with was not found for run_md5: %s" % f.get_run_md5()
        data = json.dumps(dict(success=False,
                            message=error_msg)
                          )
        return HttpResponse(data, content_type='application/json')

    # --------------------------------------
    # Was an update already sent?
    # --------------------------------------
    if dataset_run.result_received:
        error_msg = "A result was already received for this run_md5: %s" % f.get_run_md5()
        data = json.dumps(dict(success=False,
                            message=error_msg)
                          )
        return HttpResponse(data, content_type='application/json', status=400)


    # --------------------------------------
    # Update and save the DatasetScriptRun
    # --------------------------------------
    dataset_run.result_received = True  # result received
    dataset_run.result_success = f.was_run_successful() # success?
    if f.get_result_data():
        dataset_run.result_data = f.get_result_data()   # result data if available
    dataset_run.save()  # save the DatasetScriptRun

    # --------------------------------------
    # Update and save the Dataset status
    # --------------------------------------
    if dataset_run.result_success:
        dataset_run.dataset.set_status_process_success()
    else:
        dataset_run.dataset.set_status_process_failed()

    data = json.dumps(dict(success=True,
                            message="Dataset run updated")
                          )

    # Send a notification email
    send_dataset_run_message_to_tb_admins(dataset_run)

    return HttpResponse(data, content_type='application/json')
'''
import requests
payload = dict(success=True,
        run_md5='afde133c98aac9657e318de2774e687e',
        result_data='not bad')
url = 'http://127.0.0.1:8000/predict/my-dataset-run-notification/'
r = requests.post(url, data=payload)
print(r.status_code)
print(r.text)

'''