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
"""URLs for mutation data"""

from django.conf.urls import include, url

from .views import DropDownData, UploadData, UploadList, UploadView, MutationView

urlpatterns = [ # pylint: disable=invalid-name
    url(r'^json/$', DropDownData.as_view(), name="json"),
    url(r'^upload/$', UploadData.as_view(), name="upload"),
    url(r'^upload/list/$', UploadList.as_view(), name="upload.list"),
    url(r'^upload/(?P<pk>\d+)/$', UploadView.as_view(), name="upload.view"),
    url(r'^parse/$', MutationView.as_view(), name="parse"),
]
