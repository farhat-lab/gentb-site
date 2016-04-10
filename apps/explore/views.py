"""
Any views for the explore app
"""

from __future__ import print_function
from django.views.generic import DetailView
from .models import ExploreDataFileInfo

class FirstExplorePage(DetailView):
    def get_object(self):
        return ExploreDataFileInfo.objects.filter(active=True).first()

