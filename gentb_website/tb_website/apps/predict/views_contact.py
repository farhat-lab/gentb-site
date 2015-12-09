from __future__ import print_function

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext
from django.contrib.auth.decorators import login_required

from apps.predict.models import PredictDataset, PredictDatasetStatus, PredictDatasetNote
from apps.shared_data.process_file_helper import get_process_file_results
from apps.utils.view_util import get_common_dict
from apps.predict.message_helper import send_new_dataset_message_to_tb_admins


#@login_required(login_url=reverse('view_login_page', kwargs={}))
def view_dataset_contact(request, dataset_md5):
    return HttpResponse('contact here')
    d = get_common_dict(request, 'My Dataset Contact', view_my_datasets=True)

    # Not logged in, show login message
    #
    if not request.user.is_authenticated():
        return render_to_response('predict/my_dataset_contact.html',
                             d,
                             context_instance=RequestContext(request))

    try:
        dataset = PredictDataset.objects.get(md5=dataset_md5)
    except PredictDataset.DoesNotExist:
        raise Http404('PredictDataset not found')

    if not request.user == dataset.user.user:
        raise Http404('This is not your dataset')

    d['dataset'] = dataset
    d['dataset_notes'] = PredictDatasetNote.objects.filter(dataset=dataset).all()
    d['tb_user'] = dataset.user

    return render_to_response('predict/my_dataset_contact.html',
                             d,
                             context_instance=RequestContext(request))
