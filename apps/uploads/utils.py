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
Allows the downloading of URLs from different schemes using python-requests.
"""

import os
from os.path import join, getsize
import uuid
import json
from hashlib import md5

from requests import Session
from requests.compat import urlparse, urlunparse

from requests_file import FileAdapter
from requests_ftp import FTPAdapter

from django.conf import settings
from django.core.files.storage import FileSystemStorage

def get_uuid():
    """Return the hex from a new uuid4 uuid"""
    return uuid.uuid4().hex

class ManagedUrl():
    """
    Manages a URL so that it can be fed back to a user's session safely.
    It takes a url that has a password and caches the details to disk before
    sending back a hex digest hostname which points to the cache on disk.
    """
    def __init__(self, url):
        self.url = urlparse(url, 'file')
        if self.url.scheme == 'url':
            # Load the url from the disk cache
            with open(self._cache_file(self.url.netloc.encode('utf8'))) as fhl:
                self.url = self.url._replace(**json.loads(fhl.read()))

    @staticmethod
    def _cache_file(digest):
        """Return the location of the cached url"""
        return os.path.join(settings.UPLOAD_CACHE_ROOT, digest + '.json')

    def __str__(self):
        """Return the private url, what we keep internally"""
        return urlunparse(self.url)

    def __getattr__(self, name):
        return getattr(self.url, name)

    def public_url(self):
        """
        Return the public url, often the same as the private url, but
        it can be replaced by a url:// wrapper which caches passwords
        onto the disk.
        """
        url = self.url
        if url.password:
            digest = md5(url.netloc).hexdigest()
            with open(self._cache_file(digest), 'w') as fhl:
                fhl.write(json.dumps({'netloc': url.netloc, 'scheme': url.scheme}))
                url = url._replace(netloc=digest, scheme='url')
        return urlunparse(url)

    def file(self, filename):
        """Gets a managed url for a sub-file portion"""
        path = os.path.join(self.url.path, filename)
        return ManagedUrl(urlunparse(self.url._replace(path=path)))

    name = property(lambda self: os.path.basename(self.path))



class Download(): 
    """Wrap the requests module for django's storage backend."""
    _requests = None

    @classmethod
    def get_requests(cls):
        if cls._requests is None:
            cls._requests = Session()
            cls._requests.mount('file://', FileAdapter())
            cls._requests.mount('ftp://', FTPAdapter())
            #cls._requests.mount('ftps://', FTPAdapter(tls=True))
        return cls._requests

    def get_io(self, url, **kw):
        if url.username:
            kw['auth'] = (url.username, url.password)
        return self.get_requests().get(str(url), stream=True, **kw)

    def __init__(self, url, **kw):
        self.filepath = None
        self.url = ManagedUrl(url)

    def __iter__(self):
        return getattr(self, self.url.scheme.rstrip('s') + '_iter')()

    def http_iter(self):
        """List contents of a http request (file index)"""
        from bs4 import BeautifulSoup
        with self.get_io(self.url) as io:
            if io.status_code == 200:
                if io.headers.get('content-type') == 'text/html':
                    page = BeautifulSoup(io.text, 'html.parser')
                    page.prettify()
                    for anchor in page.find_all('a', href=True):
                        ret = self.url.file(anchor['href'])
                        ret.size = -1 # Don't ask for the size right now
                        yield ret
                else:
                    self.url.size = io.headers['Content-length']
                    yield self.url
            else:
                raise ValueError("Bad http response. ERR:%d" % io.status_code)

    def file_iter(self):
        """List contents of a directory"""
        path = self.url.path
        if os.path.isdir(path):
            for fn in os.listdir(path):
                full = os.path.join(path, fn) 
                if os.path.isfile(full):
                    ret = self.url.file(fn)
                    ret.size = os.path.getsize(full)
                    yield ret
        elif os.path.isfile(path):
            self.url.size = os.path.getsize(path)
            yield self.url

    def is_ok(self):
        return self.get_io(self.url).status_code == 200

    def get_error(self):
        return self.get_io(self.url).text

    def save(self, path, filename):
        """Perform the download in chunks"""
        storage = FileSystemStorage(location=path)
        self.filepath = join(path, storage.save(filename, self))

    @property
    def size(self):
        if not self.filepath:
            return 0
        return getsize(self.filepath)

    def chunks(self):
        """Translate django storage backend call, into Response iter"""
        return self.get_io(self.url).iter_content(50 * 1024) # 50KB chunks


