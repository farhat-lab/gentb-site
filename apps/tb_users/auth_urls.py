"""

 The tb_users/urls.py module uses a url prefix 'users', so login becomes
 'users:login' in a reverse match. The django 1.8 auth views do not use or
 allow one to use this prefix, so these urls will be included into
 tb_website/urls.py seperately.

"""

from django.conf.urls import patterns, url, include

def url_tree(regex, view='', *urls):
    return url(regex, include(patterns(view, *urls)))

UIDB = r'^(?P<uidb64>[0-9A-Za-z_\-]+?)/(?P<token>.+)/$'

urlpatterns = patterns('django.contrib.auth.views',
  url_tree(r'^password_reset/', 'django.contrib.auth.views',
    url(r'^$',      'password_reset',          name='password_reset'),
    url(UIDB,       'password_reset_confirm',  name='password_reset_confirm'), 
    url(r'^done/$', 'password_reset_complete', name='password_reset_complete'),
    url(r'^sent/$', 'password_reset_done',     name='password_reset_done'),
  )
)

