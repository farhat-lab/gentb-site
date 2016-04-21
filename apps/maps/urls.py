from django.conf.urls import patterns, include, url

from .views import MapPage

urlpatterns = patterns('',
    url(r'^$', MapPage.as_view(), name="map"),
)
