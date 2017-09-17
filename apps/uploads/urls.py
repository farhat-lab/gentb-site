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
from django.conf.urls import patterns, url, include 
from .views import ResumableUploadView, RetryUpload, ManualUploadView

def url_tree(regex, *urls):
    """Quick access to stitching url patterns"""
    return url(regex, include(patterns('', *urls)))

urlpatterns = patterns('', 
    url(r'^resumable/$', ResumableUploadView.as_view(), name='resumable'),
    url(r'^manual/$', ManualUploadView.as_view(), name='manual'),
    url_tree(r'^(?P<pk>\d+)/',
      url(r'^retry/$', RetryUpload.as_view(), name='retry'),
    ),
)

