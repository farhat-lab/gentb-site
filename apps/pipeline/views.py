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
Provides the views for testing and reviewing pipelines.
"""

from django.views.generic import DetailView

from .method.fake import JobManager as FakeManager
from .models import Pipeline

class PipelineDetail(DetailView):
    model = Pipeline

    def get_context_data(self, **kw):
        data = super(PipelineDetail, self).get_context_data(**kw)
        obj = self.get_object()
        data['run'] = self.get_object().run('pipeline_run',
            job_manager=FakeManager(), output_dir='~/',
            for_test=True, commit=False, file='file')
        return data

