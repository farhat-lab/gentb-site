"""
Predict app's urls
"""
#
# pylint: disable=bad-whitespace
#
from django.conf.urls import patterns, include, url

from .views import *

def url_tree(regex, *urls):
    """Quick access to stitching url patterns"""
    return url(regex, include(patterns('', *urls)))

urlpatterns = patterns('',
  url(r'^$', Datasets.as_view(), name="view_my_datasets"),
  url_tree(r'^upload/',
    url(r'^$', UploadChoices.as_view(), name="upload"),
    url(r'^(?P<type>[\w-]+)/$', UploadView.as_view(), name="upload"),
  ),
  url_tree(r'^(?P<slug>\w{32})/',
    url(r'^$',          DatasetView.as_view(),   name="view_single_dataset"),
    url(r'^note/$',     AddNote.as_view(),       name="add_note"),
  ),
)

