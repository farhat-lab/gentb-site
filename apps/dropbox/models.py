
import os
import logging
logger = logging.getLogger('apps.dropbox.models')

from urlparse import urlparse

from django.db.models import *
from django.utils.timezone import now

from .utils import Download


class DropboxFile(Model):
    url = URLField()
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

    def __str__(self):
        return str(self.filename)

    @property
    def fullpath(self):
        return os.path.join(self.file_directory, self.filename)

    def size_done(self):
        if os.path.isfile(self.fullpath):
            return os.path.getsize(self.fullpath)
        return 0

    def percent_done(self):
        """Returns the percent done of the download"""
        return int(self.size_done() / float(self.size) * 100)

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
            download.save(self.file_directory, self.filename)
            self.retrieval_end = now()
            self.size = download.size
        except Exception as error:
            self.retrieval_error = str(error)

        self.save()

