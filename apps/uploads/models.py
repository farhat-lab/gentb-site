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

from django.conf import settings
from django.db.models import (
    Model, SlugField, CharField, PositiveIntegerField, URLField,
    DateTimeField, TextField, ForeignKey
)
from django.utils.timezone import now

from .utils import Download, get_uuid
from .files import ResumableFile

LOGGER = logging.getLogger('apps.uploads.models')

class UploadFile(Model):
    """Base upload file, when a user uploads a genetic file."""
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

    # Customisable flag for apps that use this uploader
    flag = CharField(max_length=4, null=True, blank=True)

    class Meta:
        """Order by created date"""
        ordering = ('-created', 'filename')

    @classmethod
    def build_upload(cls, prefix, datum):
        """Make a new upload based on the datum"""
        obj = cls(
            name=prefix,
            filename=datum['name'],
            size=datum['bytes'],
            icon=datum['icon'],
        )
        if 'link' in datum:
            obj.url = datum['link']
        return obj

    def conclude_upload(self, directory, user=None): # pylint: disable=unused-argument
        """Finish a basic upload, should be called from child"""
        self.file_directory = directory
        return self.save()

    def get_type(self):
        """Returns the super type of any upload if possible"""
        for uname, uploader in UPLOADERS.items():
            if uname and uploader.objects.filter(pk=self.pk):
                return uname
        return 'old'

    def __str__(self):
        return "{} ({})".format(self.filename, self.get_type())

    @property
    def fullpath(self):
        """Returns the full path on disk of this file"""
        return os.path.join(self.file_directory, self.filename)

    @property
    def is_file(self):
        """Returns True if this file exists"""
        return os.path.isfile(self.fullpath)

    def size_done(self):
        """Returns the size of the file we have on disk in bytes"""
        if self.is_file:
            return os.path.getsize(self.fullpath)
        return 0

    def count(self):
        """Return a count of the number of lines in the file"""
        with open(self.fullpath, 'r') as fhl:
            return len(list(fhl))

    def percent_done(self):
        """Returns the percent done of the download"""
        return int(self.size_done() / float(self.size) * 100)

    def delete_now(self):
        """Remove the saved file from the disk if possible"""
        if self.filename and self.is_file:
            os.unlink(self.fullpath)

    def save_now(self, data):
        """Save the data as if this external download was done"""
        if not os.path.exists(self.file_directory):
            os.makedirs(self.file_directory)

        with open(self.fullpath, 'wb') as fhl:
            fhl.write(data)

        self.retrieval_start = now()
        self.retrieval_end = now()
        self.size = len(data)
        self.save()

    def download_now(self):
        """
        Download the dropbox link offline.
        """
        if not os.path.exists(self.file_directory):
            os.makedirs(self.file_directory)

        self.retrieval_start = now()
        self.retrieval_error = ''
        self.save()

        if not hasattr(self, 'url'):
            raise AttributeError("URL is a required value for downloading.")

        if self.filename is None:
            self.filename = self.url.split('/')[-1]

        try:
            download = Download(self.url)
            if download.is_ok():
                download.save(self.file_directory, self.filename)
                self.retrieval_end = now()
                self.size = download.size
            else:
                self.retrieval_error = download.get_error()

        except Exception as error: # pylint: disable=broad-except
            self.retrieval_error = str(error)
        self.save()


class DropboxUploadFile(UploadFile):
    """File uploaded via Dropbox"""
    url = URLField()


class ResumableUploadFile(UploadFile):
    """An upload from a desktop computer"""
    upload_id = SlugField(default=get_uuid)
    user = ForeignKey(settings.AUTH_USER_MODEL)

    def resumable_file(self, **kwargs):
        """Returns the resumable file object for this uf"""
        kwargs.update({
            'resumableTotalSize': self.size,
            'resumableFilename': self.filename,
        })
        return ResumableFile(self.user, kwargs)

    def conclude_upload(self, directory, user=None):
        """Actually finish off this uploaded file"""
        self.user = user
        super(ResumableUploadFile, self).conclude_upload(directory)
        self.save_resumable()

    def save_resumable(self):
        """Collect the resumable pieces and conclude it"""
        rfile = self.resumable_file()
        if rfile.is_complete:
            self.retrieval_start = rfile.started
            self.retrieval_end = rfile.ended
            self.retrieval_error = ''
            rfile.save_to(self.file_directory)
        else:
            self.retrieval_error = "File didn't completely upload"
        self.save()


class ManualUploadFile(UploadFile):
    """Uploads for sftp, local link or http"""
    url = URLField()



UPLOADERS = dict([
    (name.replace('UploadFile', '').lower(), val)
    for (name, val) in locals().items()
    if inspect.isclass(val) and issubclass(val, UploadFile)
])
