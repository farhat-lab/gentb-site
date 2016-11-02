"""
Login and user views which use the 'users:' prefix.
"""

from django.conf.urls import patterns, url

from .forms import LoginForm
from .views import UserSignUp, SignUpSuccess

urlpatterns = patterns('django.contrib.auth.views',
  url(r'^signup/$', UserSignUp.as_view(), name="signup"),
  url(r'^signup/success/$', SignUpSuccess.as_view(), name="signup_success"),

  url(r'^login/', 'login', {'authentication_form': LoginForm}, name='login'),
  url(r'^logout/', 'logout', name='logout'),
)

