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
Provide some basic front end views for pipeline testing and review.
"""

from django.urls import path

from .views import (
    PipelineDetail, PipelineRunList, PipelineRunDetail,
    JobViewer, ProgramRunDetail, ProgramRunReTry, PipelineList,
)

app_name = 'pipeline'
urlpatterns = [
    path('', PipelineList.as_view(), name='pipelines'),
    path('<int:pk>/', PipelineDetail.as_view(), name='detail'),
    path('<int:pk>/run/', PipelineRunList.as_view(), name='runs'),
    path('run/', PipelineRunList.as_view(), name='runs'),
    path('run/<int:pk>/', PipelineRunDetail.as_view(), name='run'),
    path('jobs/', JobViewer.as_view(), name='jobs'),
    path('jobs/<int:pk>/', ProgramRunDetail.as_view(), name='job'),
    path('jobs/<int:pk>/retry/', ProgramRunReTry.as_view(), name='job.retry'),
]
