from django.conf.urls import patterns, include, url

from .views import Heatmap

urlpatterns = patterns('apps.predict.views',

    url(r'^my-datasets/$', 'view_my_datasets', name="view_my_datasets"),

    url(r'^my-dataset-detail/(?P<dataset_md5>\w{32})/$', 'view_single_dataset', name="view_single_dataset"),

    url(r'^my-dataset-prediction/(?P<dataset_md5>\w{32})/$', Heatmap.as_view(), name="view_prediction_heatmap"),

)

urlpatterns += patterns('apps.predict.views_upload',

    url(r'^upload-data/$', 'view_predict_page', name="view_predict_page"),

    url(r'^upload-data-step2/(?P<dataset_md5>\w{32})/$', 'view_predict_upload_step2_confirm',
    name="view_predict_upload_step2_confirm"),

    url(r'^upload-delete/$', 'view_predict_upload_delete', name="view_predict_upload_delete"),

    url(r'^upload-success/(?P<dataset_md5>\w{32})/$', 'view_predict_upload_success', name="view_predict_upload_success"),

    #url(r'^test-upload-success/$', 'view_test_upload_success', name="view_test_upload_success"),

)

urlpatterns += patterns('apps.predict.views_run_script',
    #url(r'^my-dataset-run-script/(?P<dataset_md5>\w{32})/$', 'view_run_dataset_script', name="view_run_dataset_script"),

    url(r'^pipeline-run-results-notice/$', 'view_dataset_run_notification', name="view_dataset_run_notification"),

)
