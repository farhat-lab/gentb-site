"""
Predict app's urls
"""
#
# pylint: disable=bad-whitespace
#
from django.conf.urls import patterns, include, url

from .views import Datasets, DatasetView, UploadView, UploadConfirm, Callback

def url_tree(regex, *urls):
    """Quick access to stitching url patterns"""
    return url(regex, include(patterns('', *urls)))

urlpatterns = patterns('',
  url(r'^$',            Datasets.as_view(),      name="view_my_datasets"),
  url(r'^upload/$',     UploadView.as_view(),    name="upload"),
  url_tree(r'^(?P<slug>\w{32})/',
    url(r'^$',          DatasetView.as_view(),   name="view_single_dataset"),
    url(r'^confirm/$',  UploadConfirm.as_view(), name="view_predict_upload_step2_confirm"),
    url(r'^callback/$', Callback.as_view(),      name="view_dataset_run_notification"),
  ),
)

