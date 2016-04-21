"""
Any and all explore urls.
"""

from django.conf.urls import patterns, include, url

from .views import FirstExplorePage

urlpatterns = patterns('',
    url(r'^$', FirstExplorePage.as_view(), name="home"),
)
