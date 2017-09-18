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

from md5 import md5

from django.core.exceptions import PermissionDenied
from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required
from django.views.generic.detail import View, SingleObjectMixin
from django.views.generic import RedirectView
from django.http import HttpResponse, JsonResponse
from django.conf import settings

from .files import ResumableFile
from .models import UploadFile, ManualUploadFile
from .utils import ftp

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

        HttpResponse(status=404)

    def get_files(self, prot, url):
        """
        Sort out what protocol we should be using.
        """
        # Put some limits on the protocol size
        prot = prot.replace('_', '')[:5]
        url_hash = ManualUploadFile.cache_url(prot, url)

        for (fn, size) in getattr(self, '_' + prot)(url):
            yield {
              'id': md5(fn).hexdigest(),
              'name': fn,
              'bytes': size,
              'link': 'url://%s/%s' % (url_hash, fn),
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

    def _file(self, path):
        """
        Try and load a local directory name, if available.
        """
        # We limit users to only those with this permission
        if not self.request.user.has_perm('uploads.add_manualuploadfile'):
            raise PermissionDenied

        if os.path.isdir(path):
            for fn in os.listdir(path):
                full = os.path.join(path, fn)
                if os.path.isfile(full) and self.match_file(fn):
                    yield (fn, os.path.getsize(full))

        elif os.path.isfile(path): # Don't filter with match_file
            yield ('.', os.path.getsize(full))

    def _ftp(self, url, tls=False):
        """Load an FTP url"""
        for (name, size) in ftp(url, tls=tls):
            if self.match_file(name):
                yield (name, size)

    def _http(self, url, tls=True):
        """Load a HTTP url"""
        raise AttributeError("Protocol not supported Yet")

    def _sftp(self, url):
        return self._ftp(url, tls=True)

    def _https(self, url):
        return self._http(url, tls=True)


