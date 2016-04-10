from __future__ import print_function

from django.views.generic import TemplateView

class MapPage(TemplateView):
    template_name = 'maps/basic_map.html'

