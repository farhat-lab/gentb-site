"""
Login and user views
"""

from django.conf.urls import patterns, url
from django.contrib.auth.views import (password_reset, password_reset_done,
        password_reset_confirm, password_reset_complete)

from .forms import LoginForm

urlpatterns = patterns('apps.tb_users.views',
    url(r'^signup/$', 'view_signup_page', name="view_signup_page"),
    url(r'^signup-success/(?P<tb_user_md5>\w{32})/$',
        'view_signup_success', name="view_signup_success"),

    url(r'^accounts/password_reset/$', password_reset,
        {'post_reset_redirect' : '/accounts/password_reset/mailed/'},
        name="password_reset"),

    (r'^accounts/password_reset/mailed/$', password_reset_done),

    url(r'^accounts/password_reset/(?P<uidb64>[0-9A-Za-z_\-]+)/'
        r'(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        password_reset_confirm,
        {'post_reset_redirect' : '/accounts/password_reset/complete/'},
        name="password_reset_confirm"),

    (r'^accounts/password_reset/complete/$', password_reset_complete),
)

urlpatterns += patterns('django.contrib.auth.views',
  url(r'^login/', 'login', {'authentication_form': LoginForm}, name='login'),
  url(r'^logout/', 'logout', name='logout'),
)

