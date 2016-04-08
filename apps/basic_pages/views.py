from __future__ import print_function

from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext
from apps.explore.models import ExploreDataFileInfo
from django.views.generic import TemplateView



def view_explore_page(request):

    d = get_common_dict(request, 'Explore', explore_page=True)

    explore_link_info = ExploreDataFileInfo.objects.filter(active=True).first()
    d.update(dict(explore_link_info=explore_link_info))

    return render_to_response('explore.html'\
                            , d\
                            , context_instance=RequestContext(request))


