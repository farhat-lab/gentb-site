
import requests
import logging
logger = logging.getLogger('apps.dropbox.models')

from StringIO import StringIO
from urlparse import urlparse

from django.utils.timezone import now
from django.core.files import File
from django.db.models import *

from apps.predict.models import PredictDataset

class DropboxFile(Model):
    dataset = ForeignKey(PredictDataset, related_name='files')

    url = URLField()
    name = SlugField(max_length=32)
    result = FileField(null=True, blank=True)
    created = DateTimeField(auto_now_add=True)

    # system attempts to download files
    retrieval_start = DateTimeField(null=True, blank=True)
    retrieval_end = DateTimeField(null=True, blank=True)
    retrieval_error = TextField(blank=True)

    class Meta:
        ordering = ('-created', 'dataset')

    def __str__(self):
        return '{0} ({1})'.format(self.dataset, self.name)

    @property
    def filename(self):
        return urlparse(self.url).path.split('/')[-1]

    def download_now(self):
        """
        Download the dropbox link offline.
        """
        self.retrieval_start = now()
        self.retrieval_error = ''
        self.save()

        try:
            io = StringIO(requests.get(self.url).content)
            self.result =  File(io, name=self.filename)
            self.retrieval_end = now()
        except Exception as error:
            self.retrieval_error = str(error)
            raise
        self.save()

