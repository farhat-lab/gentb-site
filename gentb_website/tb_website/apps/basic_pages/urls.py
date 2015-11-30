from django.conf.urls import patterns, include, url


urlpatterns = patterns('apps.basic_pages.views',

    url(r'^home/$', 'view_homepage', name="view_homepage"),


    url(r'^about/$', 'view_about_page', name="view_about_page"),

    #url(r'^predict/$', 'view_predict_page', name="view_predict_page"),

    url(r'^share/$', 'view_data_page', name="view_data_page"),

    url(r'^data-upload-information/$', 'view_data_upload_information', name="view_data_upload_information"),

    url(r'^explore/$', 'view_explore_page', name="view_explore_page"),

    url(r'^/?$', 'view_homepage', name="default_homepage"),

    #url(r'^milestone-history/(?P<chosen_year>(\d){4})/$', 'view_milestone_history', name="view_milestone_history_by_year"),

    #url(r'^milestone-roadmap/(?P<repo_name>(\-|_|\w){1,120})/$', 'view_single_repo_column', name="view_single_repo_column"),

)
