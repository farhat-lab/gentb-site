from django.conf.urls import patterns, include, url

from .views import *

urlpatterns = patterns('',
    url(r'^$', MapPage.as_view(), name="map"),
    url(r'data/places/$',    Places.as_view(), name="map.places"),
    url(r'data/drugs/$',     Drugs.as_view(), name="map.drugs"),
    url(r'data/lineages/$',  Lineages.as_view(), name="map.lineages"),
    url(r'data/mutation/$',  MutationView.as_view(), name="map.mutation"),
    url(r'data/mutations/$', Mutations.as_view(), name="map.mutations"),
)
