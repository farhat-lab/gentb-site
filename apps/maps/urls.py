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
# pylint: disable=invalid-name
"""
Maps urls
"""

from django.urls import path
from django.views.generic.base import RedirectView

from .views import MapPage, Places, DrugList, Lineages, \
    MutationView, LocusList, Mutations, Sources

from .views_antibiograms import AntibiogramMap, MarginalPlaces, MarginalDrugs

app_name = 'maps'
urlpatterns = [
    # Default map, used for quick url access
    path('', RedirectView.as_view(pattern_name='map.mutations'), name='map'),

    # Original mutations map
    path('mutations/', MapPage.as_view(), name="map.mutations"),
    path('data/sources/', Sources.as_view(), name="data.sources"),
    path('data/places/', Places.as_view(), name="data.places"),
    path('data/drugs/', DrugList.as_view(), name="data.drugs"),
    path('data/lineages/', Lineages.as_view(), name="data.lineages"),
    path('data/locuses/', LocusList.as_view(), name="data.locuses"),
    path('data/mutations/', Mutations.as_view(), name="data.mutations"),
    path('data/mutation/', MutationView.as_view(), name="data.mutation"),

    # New autobiograms map
    path('antibiogram/', AntibiogramMap.as_view(), name="map.antibiogram"),
    path('data/antibiograms/places/', MarginalPlaces.as_view(), name="data.marginalplaces"),
    path('data/antibiograms/drugs/', MarginalDrugs.as_view(), name="data.marginaldrugs"),

]
