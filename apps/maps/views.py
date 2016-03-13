from __future__ import print_function

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext

from apps.utils.view_util import get_common_dict


def view_map_page(request):

  d = get_common_dict(request, 'Map', map_page=True)

  return render_to_response('maps/basic_map.html'\
                            , d\
                            , context_instance=RequestContext(request))
