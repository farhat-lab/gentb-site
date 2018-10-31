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
Provide some basic front end views for pipeline testing and review.
"""

from django.conf.urls import url

from .views import (
    PipelineDetail, PipelineRunList, PipelineRunDetail,
    JobViewer, ProgramRunDetail, ProgramRunReTry, PipelineList,
)

urlpatterns = [ # pylint: disable=invalid-name
    url(r'^$', PipelineList.as_view(), name='pipelines'),
    url(r'^(?P<pk>\d+)/$', PipelineDetail.as_view(), name='detail'),
    url(r'^(?P<pk>\d+)/run/$', PipelineRunList.as_view(), name='runs'),
    url(r'^run/$', PipelineRunList.as_view(), name='runs'),
    url(r'^run/(?P<pk>\d+)/$', PipelineRunDetail.as_view(), name='run'),
    url(r'^jobs/$', JobViewer.as_view(), name='jobs'),
    url(r'^jobs/(?P<pk>\d+)/$', ProgramRunDetail.as_view(), name='job'),
    url(r'^jobs/(?P<pk>\d+)/retry/$', ProgramRunReTry.as_view(), name='job.retry'),
]
