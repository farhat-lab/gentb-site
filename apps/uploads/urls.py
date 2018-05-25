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
"""
Urls to upload files to the system.
"""
from django.conf.urls import url
from .views import ResumableUploadView, RetryResumableUpload, \
    RetryUpload, ManualUploadView, TestUpload

urlpatterns = [
    url(r'^test/$', TestUpload.as_view(), name='test'),
    url(r'^resumable/$', ResumableUploadView.as_view(), name='resumable'),
    url(r'^resumable/(?P<pk>\d+)/$', RetryResumableUpload.as_view(), name='resumable'),
    url(r'^manual/$', ManualUploadView.as_view(), name='manual'),
    url(r'^(?P<pk>\d+)/retry/$', RetryUpload.as_view(), name='retry'),
]
