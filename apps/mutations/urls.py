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
URLs for mutation data
"""

from django.urls import path

from .views import DropDownData, UploadData, UploadList, UploadView, MutationView

app_name = 'genes'
urlpatterns = [ # pylint: disable=invalid-name
    path('json/', DropDownData.as_view(), name="json"),
    path('upload/', UploadData.as_view(), name="upload"),
    path('upload/list/', UploadList.as_view(), name="upload.list"),
    path('upload/<int:pk>/', UploadView.as_view(), name="upload.view"),
    path('parse/', MutationView.as_view(), name="parse"),
]
