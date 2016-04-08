from django.conf.urls import patterns, include, url


urlpatterns = patterns('apps.tb_users.views_signup',
    url(r'^signup/$', 'view_signup_page', name="view_signup_page"),
    url(r'^signup-success/(?P<tb_user_md5>\w{32})/$', 'view_signup_success', name="view_signup_success"),

)

urlpatterns += patterns('apps.tb_users.views_login',
    url(r'^login/$', 'view_login_page', name="view_login_page"),
    url(r'^logout/$', 'view_logout_page', name="view_logout_page"),
)
