from __future__ import print_function

from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext
from apps.utils.view_util import get_common_dict
from apps.explore.models import ExploreDataFileInfo


def view_homepage(request, just_logged_in=False):

    d = get_common_dict(request, 'Home', home_page=True)

    d['just_logged_in'] = just_logged_in

    #d['page_title'] = 'Phthisis Ravens: TB Project'

    return render_to_response('homepage.html'\
                              , d\
                              , context_instance=RequestContext(request))


def view_about_page(request):
    d = get_common_dict(request, 'About', about_page=True)

    return render_to_response('about.html'\
                              , d\
                              , context_instance=RequestContext(request))


def view_share_page(request):
    d = get_common_dict(request, 'Share', share_page=True)

    return render_to_response('share.html'\
                            , d\
                            , context_instance=RequestContext(request))




def view_explore_page(request):

    d = get_common_dict(request, 'Explore', explore_page=True)

    explore_link_info = ExploreDataFileInfo.objects.filter(active=True).first()
    d.update(dict(explore_link_info=explore_link_info))

    return render_to_response('explore.html'\
                            , d\
                            , context_instance=RequestContext(request))
