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

from django.views.generic import ListView, DetailView, TemplateView

from chore.fake import FakeJobManager
from chore import get_job_manager
from .models import Pipeline, PipelineRun, ProgramRun

class PipelineDetail(DetailView): # pylint: disable=too-many-ancestors
    """Test the pipeline in the front end to test the connectivity"""
    model = Pipeline

    def get_context_data(self, **kw):
        data = super(PipelineDetail, self).get_context_data(**kw)
        self.get_object()
        data['run'] = self.get_object().run('pipeline_run',\
            job_manager=FakeJobManager(), output_dir='~/',\
            for_test=True, commit=False, file='file')
        return data

class PipelineRunList(ListView): # pylint: disable=too-many-ancestors
    """Create a list of run pipelines and their results"""
    model = PipelineRun
    paginate_by = 10
    ordering = ['-created']

class PipelineRunDetail(DetailView): # pylint: disable=too-many-ancestors
    """Show a single pipeline run and it's component runs"""
    model = PipelineRun

class ProgramRunDetail(DetailView): # pylint: disable=too-many-ancestors
    """Show a single program run"""
    model = ProgramRun

class JobViewer(TemplateView):
    """Lets users view a list of jobs and how they are running"""
    template_name = 'pipeline/jobs_status_list.html'

    def get_context_data(self, **kw):
        data = super(JobViewer, self).get_context_data(**kw)
        data['pipeline'] = get_job_manager()
        kw = {}
        if 'user' in self.request.GET:
            kw['user'] = self.request.GET['user']
        cols = [c for c in self.request.GET.get('cols', '').split(',') if c]
        data['object_list'] = [self.get_item(item, cols)\
            for item in data['pipeline'].jobs_status(*cols, **kw)]
        data['pipeline_name'] = type(data['pipeline']).__name__
        data['cols'] = cols
        return data

    @staticmethod
    def get_item(item, cols):
        """Parses the data from item for view"""
        finished = 'success'
        if item['error'] and item['error'] != 'None':
            finished = 'danger'
        item['cols'] = [item.get(c, None) for c in cols]
        item['context'] = {
            'finished': finished,
            'pending': 'default',
            'running': 'primary',
        }.get(item['status'].lower(), 'default')
        return item
