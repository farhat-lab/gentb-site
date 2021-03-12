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
Control the resumable file API (javascript client side)
"""

import os
from datetime import datetime

import pytz

from django.core.exceptions import ImproperlyConfigured
from django.core.files.storage import FileSystemStorage
from django.conf import settings

class ResumableFile(object):
    """A resumable file controls getting pieces of a file"""
    upload_root = getattr(settings, 'UPLOAD_ROOT', None)
    name_template = "part_%(resumableChunkNumber)04d.tmp"

    def __init__(self, user, kwargs):
        self.kwargs = kwargs
        self.user = user

        if not self.upload_root:
            raise ImproperlyConfigured('You must set UPLOAD_ROOT in settings')

        self.storage = FileSystemStorage(location=self.upload_dir)

    @property
    def chunk_exists(self):
        """Checks if the requested chunk exists."""
        name = self.name_template % self.kwargs
        if not self.storage.exists(name):
            return False
        chunk_size = int(self.kwargs.get('resumableCurrentChunkSize'))
        return self.storage.size(name) == chunk_size

    def chunk_names(self):
        """Iterates over all stored chunks and yields their names."""
        try:
            return sorted(self.storage.listdir('')[1])
        except OSError:
            return []

    def chunks(self):
        """Yield the contents of every chunk, FileSystemStorage.save compatible"""
        for name in self.chunk_names():
            try:
                yield self.storage.open(name).read()
            except AttributeError:
                raise IOError("Couldn't read {}".format(name))

    def delete_chunks(self):
        """Remove every chunk (once complete)"""
        return [self.storage.delete(chunk) for chunk in self.chunk_names()]

    @property
    def filename(self):
        """Gets the filename."""
        filename = self.kwargs.get('resumableFilename')
        if '/' in filename:
            raise Exception('Invalid filename')
        return filename

    @property
    def upload_dir(self):
        """Gets the directory to save chunks to"""
        return os.path.join(self.upload_root, str(self.user.pk), self.filename)

    @property
    def is_complete(self):
        """Checks if all chunks are allready stored."""
        return os.path.isfile(self.upload_dir) or \
            self.kwargs['resumableTotalSize'] == self.size

    def process_chunk(self, _file):
        """Process the chunks of the given file"""
        if not self.chunk_exists:
            try:
                self.storage.save(self.name_template % self.kwargs, _file)
            except AttributeError:
                pass # Error saving file
            #except (IOError, OSError) as err:
                #raise IOError("Tried to save: {}, {}, {}".format(self.name_template, self.kwargs, self.name_template % self.kwargs))
                #pass # Existing file in the way

    def save_to(self, new_dir):
        """When saving all the chunks to a new directory"""
        filename = os.path.join(new_dir, self.filename)

        if os.path.islink(self.upload_dir):
            # This was previously uploaded and we can relink it.
            if not os.path.exists(filename):
                linkto = os.readlink(self.upload_dir)
                os.symlink(linkto, filename)
            return

        # Actually save the file using storage
        storage = FileSystemStorage(location=new_dir)
        storage.save(self.filename, self)

        # Delete all the chunks after use
        self.delete_chunks()

        # Sometimes there's arace condition if people double-click
        if os.path.isdir(self.upload_dir) and not os.listdir(self.upload_dir):
            os.rmdir(self.upload_dir)

        # Create a symlink for tracking and re-user
        if os.path.isfile(filename):
            os.symlink(filename, self.upload_dir)

    @property
    def size(self):
        """Gets chunks size."""
        size = 0
        for chunk in self.chunk_names():
            size += self.storage.size(chunk)
        return size

    @property
    def started(self):
        """Return the first modified datetime"""
        return self.get_times()[0]

    @property
    def ended(self):
        """Return the last modified datetime"""
        return self.get_times()[-1]

    def get_times(self):
        """Return a list of modified datetimes"""
        upload = self.upload_dir
        files = []
        if os.path.isdir(upload):
            files = [os.path.join(upload, f) for f in os.listdir(upload)]
        elif os.path.exists(upload):
            files = [upload]
        return [fromtimestamp(os.path.getmtime(f)) for f in files]

def fromtimestamp(dtm):
    """Format a string timestamp into a non-native UTC timestamp"""
    return pytz.utc.localize(datetime.fromtimestamp(dtm))
