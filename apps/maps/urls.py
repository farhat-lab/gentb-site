#
# Copyright (C) 2017 Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from django.conf.urls import include, url

from .views import MapPage, Places, DrugList, Lineages, LineageBreakdown, \
    MutationView, LocusRange, LocusList, Mutations, Sources

urlpatterns = [
    url(r'^$', MapPage.as_view(), name="map"),
    url(r'data/places/$', Places.as_view(), name="map.places"),
    url(r'data/drugs/$', DrugList.as_view(), name="map.drugs"),
    url(r'data/lineages/$', Lineages.as_view(), name="map.lineages"),
    url(r'data/lineages/breakdown/$', LineageBreakdown.as_view(), name="map.lineage_breakdown"),
    url(r'data/locrange/$', LocusRange.as_view(), name="map.locus.range"),
    url(r'data/locuses/$', LocusList.as_view(), name="map.locuses"),
    url(r'data/mutation/$', MutationView.as_view(), name="map.mutation"),
    url(r'data/mutations/$', Mutations.as_view(), name="map.mutations"),
    url(r'data/sources/$', Sources.as_view(), name="map.sources"),
]
