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
Urls to upload files to the system.
"""
from django.urls import path

from .views import ResumableUploadView, RetryResumableUpload, \
    RetryUpload, ManualUploadView, TestUpload

app_name = 'uploads'
urlpatterns = [
    path('test/', TestUpload.as_view(), name='test'),
    path('resumable/', ResumableUploadView.as_view(), name='resumable'),
    path('resumable/<int:pk>/', RetryResumableUpload.as_view(), name='resumable'),
    path('manual/', ManualUploadView.as_view(), name='manual'),
    path('<int:pk>/retry/', RetryUpload.as_view(), name='retry'),
]
