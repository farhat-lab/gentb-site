from django.conf.urls import patterns, include, url


urlpatterns = patterns('apps.shared_data.views',

    url(r'^upload-data/$', 'view_share_page', name="view_share_page"),

    url(r'^upload-success/(?P<shared_info_md5>\w{32})/$', 'view_upload_success', name="view_upload_success"),


    #url(r'^milestone-history/(?P<chosen_year>(\d){4})/$', 'view_milestone_history', name="view_milestone_history_by_year"),

    #url(r'^milestone-roadmap/(?P<repo_name>(\-|_|\w){1,120})/$', 'view_single_repo_column', name="view_single_repo_column"),

)


