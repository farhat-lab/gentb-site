from django.conf.urls import patterns, include, url


urlpatterns = patterns('apps.tb_users.views_signup',
    url(r'^signup/$', 'view_signup_page', name="view_signup_page"),
    url(r'^signup-success/(?P<tb_user_md5>\w{32})/$', 'view_signup_success', name="view_signup_success"),

    url(r'^accounts/password_reset/$',
        'django.contrib.auth.views.password_reset',
        {'post_reset_redirect' : '/accounts/password_reset/mailed/'},
        name="password_reset"),

    (r'^accounts/password_reset/mailed/$',
        'django.contrib.auth.views.password_reset_done'),

    url(r'^accounts/password_reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        'django.contrib.auth.views.password_reset_confirm',
        {'post_reset_redirect' : '/accounts/password_reset/complete/'},
                name="password_reset_confirm"),

    (r'^accounts/password_reset/complete/$',
        'django.contrib.auth.views.password_reset_complete'),
)

urlpatterns += patterns('apps.tb_users.views_login',
    url(r'^login/$', 'view_login_page', name="view_login_page"),
    url(r'^logout/$', 'view_logout_page', name="view_logout_page"),
)
