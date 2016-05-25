#
# pylint: disable=bad-whitespace
#
from django.conf.urls import patterns, include, url

from .views import *

def url_tree(regex, *urls):
    """Quick access to stitching url patterns"""
    return url(regex, include(patterns('', *urls)))

urlpatterns = patterns('',
  # None yet
)

