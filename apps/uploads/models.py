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
Manage uploads to django via a number of different mechanisms.
"""

import os
import inspect
import logging
logger = logging.getLogger('apps.uploads.models')

from urlparse import urlparse

from django.db.models import *
from django.utils.timezone import now

from .utils import Download


class UploadFile(Model):
    name = SlugField(max_length=128)
    filename = CharField(max_length=255, null=True)
    file_directory = CharField(max_length=255)

    size = PositiveIntegerField(null=True)
    icon = URLField(null=True)

    created = DateTimeField(auto_now_add=True)

    # system attempts to download files
    retrieval_start = DateTimeField(null=True, blank=True)
    retrieval_end = DateTimeField(null=True, blank=True)
    retrieval_error = TextField(blank=True)

    class Meta:
        ordering = ('-created', 'filename')

    @classmethod
    def build_upload(cls, prefix, datum):
        return cls(
	    name=prefix,
	    filename=datum['name'],
	    size=datum['bytes'],
	    icon=datum['icon'],
        )

    def __str__(self):
        return str(self.filename)

    @property
    def fullpath(self):
        return os.path.join(self.file_directory, self.filename)

    @property
    def is_file(self):
        return os.path.isfile(self.fullpath)

    def size_done(self):
        if self.is_file:
            return os.path.getsize(self.fullpath)
        return 0

    def percent_done(self):
        """Returns the percent done of the download"""
        return int(self.size_done() / float(self.size) * 100)

    def delete_now(self):
        """Remove the saved file from the disk if possible"""
        if self.filename and self.is_file:
            os.unlink(self.fullpath)

    def save_now(self, data):
        """Save the data as if this dropbox download was done"""
        if not os.path.exists(self.file_directory):
            os.makedirs(self.file_directory)

        with open(self.fullpath, 'wb') as fhl:
            fhl.write(data)

        self.retrieval_start = now()
        self.retrieval_end = now()
        self.size = len(data)
        self.save()


class DropboxUploadFile(UploadFile):
    """File uploaded via Dropbox"""
    url = URLField()

    @classmethod
    def build_upload(cls, datum):
        obj = super(cls, cls).build_upload(datum)
	obj.url = datum['link']
        return obj

    def download_now(self):
        """
        Download the dropbox link offline.
        """
        if not os.path.exists(self.file_directory):
            os.makedirs(self.file_directory)

        self.retrieval_start = now()
        self.retrieval_error = ''
        self.save()

        if self.filename is None:
            self.filename = self.url.split('/')[-1]

        try:
            download = Download(self.url)
            if download.is_ok():
                download.save(self.file_directory, self.filename)
                self.retrieval_end = now()
                self.size = download.size
            else:
                self.retrieval_error = download.io.text

        except Exception as error:
            self.retrieval_error = str(error)

        self.save()

UPLOADERS = dict([
    (name.replace('UploadFile', '').lower(), cls)
        for (name, cls) in locals().items()
            if inspect.isclass(cls) and issubclass(cls, UploadFile)
    ])

