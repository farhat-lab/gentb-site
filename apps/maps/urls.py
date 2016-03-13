from django.conf.urls import patterns, include, url


urlpatterns = patterns('apps.maps.views',

    url(r'^tb-map/$', 'view_map_page', name="view_map_page"),

    #url(r'^(?P<shared_info_md5>\w{32})/$', 'view_upload_success', name="view_upload_success"),


)


