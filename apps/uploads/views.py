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

from md5 import md5

from django.core.exceptions import PermissionDenied
from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required
from django.views.generic.detail import View, SingleObjectMixin
from django.views.generic import RedirectView, FormView
from django.http import HttpResponse, JsonResponse

from .files import ResumableFile
from .models import UploadFile, ResumableUploadFile
from .utils import Download

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
        """Checks if chunk has allready been sent."""
        resumable = self.get_object()
        if not (resumable.chunk_exists or resumable.is_complete):
            return HttpResponse('Chunk not found', status=204)
        return HttpResponse('Chunk already exists')

    def post(self, *args, **kwargs):
        """Saves chunks then checks if the file is complete."""
        resumable = self.get_object()
        if resumable.chunk_exists:
            return HttpResponse('chunk already exists')
        resumable.process_chunk(self.request.FILES.get('file'))
        return HttpResponse()

class RetryResumableUpload(ResumableUploadView):
    def get_object(self):
        if not hasattr(self, 'upload'):
            self.upload = ResumableUploadFile.objects.get(
                pk=self.kwargs['pk'], user=self.request.user)
        return self.upload.resumable_file(**dict(self._get_kwargs()))

    def dispatch(self, request, *args, **kw):
        """Runs the chunk upload and then checks for completeness"""
        ret = super(RetryResumableUpload, self).dispatch(request, *args, **kw)
        if self.get_object().is_complete:
            self.upload.save_resumable()
        return ret

class ManualUploadView(View):
    """
    This view takes a URL and attempts to access the files that might be inside.
    """
    def post(self, *args, **kwargs):
        """Get the URL and process"""
        url = unicode(self.request.POST.get('url', ''))

        # Support local /server/ filenames
        if url.startswith('/'):
            url = 'file://' + url

        if url.startswith('file://'):
            if not self.request.user.has_perm('uploads.add_manualuploadfile'):
                raise PermissionDenied

        if '://' in url:
            try:
                if not self.request.user.is_authenticated():
                    raise PermissionDenied("Not logged in")
                return JsonResponse({
                  'files': list(self.get_files(*url.split('://', 1))),
                })
            except PermissionDenied as err:
                return HttpResponse(str(err), status=403)
            except AttributeError as err:
                return HttpResponse(str(err), status=404)

        return HttpResponse(status=404)

    def get_files(self, prot, url):
        """
        Sort out what protocol we should be using.
        """
        for url in Download(url):
            if self.match_file(url.name):
                yield {
                  'id': md5(str(url)).hexdigest(),
                  'name': url.name,
                  'bytes': url.size,
                  'link': url.public_url(),
                  # XXX To replace with better icons
                  'icon': 'https://www.dropbox.com/static/images/icons64/page_white_compressed.png',
                }

    def match_file(self, filename):
        """
        Match a filename with any filter instructions
        return True if it matches the filters.
        """
        if 'extensions' in self.request.POST:
            for ext in self.request.POST['extensions'].split(" "):
                if filename.endswith(ext):
                    return True
            return False
        return True



from .forms import TestUploadForm
class TestUpload(FormView):
    form_class = TestUploadForm
    template_name = 'uploads/test_form.html'

    def form_valid(self, form):
        """Because we are testing, don't redirect"""
        return self.render_to_response({'form': form})

