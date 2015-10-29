from django.conf.urls import patterns, include, url

urlpatterns = patterns('apps.predict.views',

    url(r'^my-datasets/$', 'view_my_datasets', name="view_my_datasets"),

    url(r'^my-dataset-detail/(?P<dataset_md5>\w{32})/$', 'view_single_dataset', name="view_single_dataset"),

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
    url(r'^my-dataset-run-script/(?P<dataset_md5>\w{32})/$', 'view_run_dataset_script', name="view_run_dataset_script"),

    url(r'^my-dataset-run-notification/$', 'view_dataset_run_notification', name="view_dataset_run_notification"),

)

urlpatterns += patterns('apps.predict.views_dropbox_download',
    url(r'^file-retrieval-results/$', 'record_file_retrieval_results', name="record_file_retrieval_results"),

)
#urlpatterns += patterns('apps.predict.views_contact',
    #url(r'^my-dataset-contact/(?P<dataset_md5>\w{32})/$', 'view_dataset_contact', name="view_dataset_contact"),
#)
