#
# pylint: disable=bad-whitespace
#
from django.conf.urls import patterns, include, url

from .views import *

urlpatterns = patterns('',
  url(r'^json/$', DropDownData.as_view(), name="json"),
  url(r'^upload/$', UploadData.as_view(), name="upload"),
  url(r'^upload/(?P<pk>\d+)/$', UploadView.as_view(), name="upload.view"),
  url(r'^parse/$', MutationView.as_view(), name="parse"),
)

