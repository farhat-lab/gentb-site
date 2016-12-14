from django.conf.urls import patterns, include, url

from .views import MapPage, Countries

urlpatterns = patterns('',
    url(r'^$', MapPage.as_view(), name="map"),
    url(r'countries/$', Countries.as_view(), name="map.countries"),
)
