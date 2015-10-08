from __future__ import print_function

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext

from apps.shared_data.forms import SharedFileInfoForm
from apps.shared_data.models import SharedFileInfo
from apps.shared_data.process_file_helper import get_process_file_results

from apps.utils.view_util import get_common_dict


def view_map_page(request):

  d = get_common_dict(request, 'Map', map_page=True)

  return render_to_response('maps/basic_map.html'\
                            , d\
                            , context_instance=RequestContext(request))

"""
def view_share_page(request):
    d = {}
    d['share_page'] = True
    
    print (request.FILES)
    if request.POST:
        f = SharedFileInfoForm(request.POST, request.FILES, label_suffix='')
        if f.is_valid():
            shared_file_info = f.save()
            success_url = reverse('view_upload_success'\
                                , kwargs=dict(shared_info_md5=shared_file_info.md5)\
                                )
            return HttpResponseRedirect(success_url)
        else:
            d['ERROR_FOUND']  = True
    else:
        f = SharedFileInfoForm(label_suffix='')
    
    d['share_form'] = f
    
    return render_to_response('share/share.html'\
                            , d\
                            , context_instance=RequestContext(request))

"""
