#
# Copyright (C) 2017 Maha Farhat
#           (C) 2015 Jean-Philippe Serafin (MIT)
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
# Imported from django-resumable by jean-philippe serafin 0.2.0-dev, (c) 2015
# MIT License (no file header or project LICENSE file, see setup.py)
#
"""
Allow uploads to be 'chunked' and saved in descrete chunks.
"""

import os
import re

from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required
from django.views.generic.detail import View, SingleObjectMixin
from django.views.generic import RedirectView
from django.http import HttpResponse
from django.conf import settings

from .files import ResumableFile
from .models import UploadFile

class RetryUpload(SingleObjectMixin, RedirectView):
    model = UploadFile
    permanent = False

    def get_redirect_url(self, **kwargs):
        upload = self.get_object()
        upload.retrieval_start = None
        upload.retrieval_end = None
        upload.retrieval_error = ''
        upload.save()
        return self.request.GET.get('next', '/')


class ResumableUploadView(View):
    @method_decorator(login_required)
    def dispatch(self, request, *args, **kw):
        return super(ResumableUploadView, self).dispatch(request, *args, **kw)

    def get_object(self):
        """Returns the ResumableFile object for this chunk"""
        return ResumableFile(self.request.user, dict(self._get_kwargs()))

    def _get_kwargs(self):
        """Returns the GET request pre-formated"""
        for name, value in self.request.GET.items():
            if value.isdigit():
                yield (name, int(value))
            elif value in ('true', 'false'):
                yield (name, value == 'true')
            else:
                yield (name, value)

    def get(self, *args, **kwargs):
        """Checks if chunk has allready been sended."""
        r = self.get_object()
        if not (r.chunk_exists or r.is_complete):
            return HttpResponse('Chunk not found', status=204)
        return HttpResponse('Chunk already exists')

    def post(self, *args, **kwargs):
        """Saves chunks then checks if the file is complete."""
        r = self.get_object()
        if r.chunk_exists:
            return HttpResponse('chunk already exists')
        r.process_chunk(self.request.FILES.get('file'))
        return HttpResponse()


class UrlUploadView(View):
    """
    This view takes a URL and attempts to access the files that might be inside.
    """
    def post(self, *args, **kwargs):
        """Get the URL and process"""
        pass

