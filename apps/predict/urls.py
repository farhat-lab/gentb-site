"""
Predict app's urls
"""
#
# pylint: disable=bad-whitespace
#
from django.conf.urls import patterns, include, url

from .views import Datasets, DatasetView, Heatmap, UploadView

def url_tree(regex, *urls):
    """Quick access to stitching url patterns"""
    return url(regex, include(patterns('', *urls)))

urlpatterns = patterns('',
  url(r'^$',           Datasets.as_view(),    name="view_my_datasets"),
  url(r'^upload/$',    UploadView.as_view(),  name="view_predict_page"),
  url_tree(r'^(?P<slug>\w{32})/',
    url(r'^$',         DatasetView.as_view(), name="view_single_dataset"),
    url(r'^heatmap/$', Heatmap.as_view(),     name="view_prediction_heatmap"),
    #url(r'^confirm/$', UploadConfirm.as_view(), name="view_predict_upload_step2_confirm"),
    #url(r'^delete/$',  DeleteUpload.as_view(), name="view_predict_upload_delete"),
  ),
)

#urlpatterns += patterns('apps.predict.views_upload',
#  DONE   url(r'^upload-data-step2/(?P<dataset_md5>\w{32})/$', 'view_predict_upload_step2_confirm', name="view_predict_upload_step2_confirm"),
#  DONE  url(r'^upload-delete/$', 'view_predict_upload_delete', name="view_predict_upload_delete"),
#  REPLACE?  url(r'^upload-success/(?P<dataset_md5>\w{32})/$', 'view_predict_upload_success', name="view_predict_upload_success"),
#)

#urlpatterns += patterns('apps.predict.views_run_script',
    # PRE_DISABLED url(r'^my-dataset-run-script/(?P<dataset_md5>\w{32})/$', 'view_run_dataset_script', name="view_run_dataset_script"),
    # UNKNOWN-FUNCTION url(r'^pipeline-run-results-notice/$', 'view_dataset_run_notification', name="view_dataset_run_notification"),
#)
