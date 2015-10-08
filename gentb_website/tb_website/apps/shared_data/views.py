from __future__ import print_function

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext

from apps.shared_data.forms import SharedFileInfoForm
from apps.shared_data.models import SharedFileInfo
from apps.shared_data.process_file_helper import get_process_file_results

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


def view_upload_success(request, shared_info_md5):

    try:
        shared_file_info = SharedFileInfo.objects.get(md5=shared_info_md5)
    except SharedFileInfo.DoesNotExist:
        raise Http404('SharedFileInfo not found')
    
    #(success, msg_or_data) = get_process_file_results(shared_file_info.file_obj.file.name)
    #print ('msg_or_data', type(msg_or_data))
    d = {}
    d['UPLOAD_SUCCESS'] = True
    
    #d['FILE_PROCESS_SUCCESS'] = success
    #d['FILE_PROCESS_ERR_OR_DATA'] = msg_or_data
    
    d['share_page'] = True
    d['shared_file_info'] = shared_file_info
    
    
    return render_to_response('share/share.html'\
                            , d\
                            , context_instance=RequestContext(request))
    
    