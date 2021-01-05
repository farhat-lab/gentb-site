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
Predict app's urls
"""

from django.urls import path
from django.conf.urls import include, url

from .views import (
    Datasets, UploadChoices, UploadView, DatasetView, AddNote, ScatterPlot,
    DatasetViewProcessing, DatasetViewOutput, DatasetViewPredict, DatasetViewLineages
)

def url_tree(regex, *urls):
    class UrlTwig():
        urlpatterns = urls
    return url(regex, include(UrlTwig))

app_name = 'predict'
urlpatterns = [
    path('', Datasets.as_view(), name="view_my_datasets"),
    url_tree(
        r'^upload/',
        url(r'^$', UploadChoices.as_view(), name="upload"),
        url(r'^(?P<type>[\w-]+)/$', UploadView.as_view(), name="upload"),
    ),
    url_tree(
        r'^(?P<slug>\w{32})/',
        url(r'^$', DatasetView.as_view(), name="view_single_dataset"),
        url_tree(
            r'^page/',
            url(r'^process/$', DatasetViewProcessing.as_view(), name="dataset_proc"),
            url(r'^output/$', DatasetViewOutput.as_view(), name="dataset_out"),
            url(r'^predict/$', DatasetViewPredict.as_view(), name="dataset_pred"),
            url(r'^lineages/$', DatasetViewLineages.as_view(), name="dataset_lin"),
        ),
        url(r'^note/$', AddNote.as_view(), name="add_note"),
    ),
    url(r'^results/(?P<pk>\d+)/plot/', ScatterPlot.as_view(), name="scatter_plot"),
]
