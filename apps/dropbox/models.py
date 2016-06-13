
import logging
logger = logging.getLogger('apps.dropbox.models')

from urlparse import urlparse

from django.db.models import *
from django.utils.timezone import now

from apps.predict.models import PredictDataset, PredictDatasetFile

from .utils import Download


class DropboxFile(Model):
    dataset = ForeignKey(PredictDataset, related_name='files')
    result = ForeignKey(PredictDatasetFile, null=True, blank=True)

    url = URLField()
    name = SlugField(max_length=32)
    created = DateTimeField(auto_now_add=True)

    # system attempts to download files
    retrieval_start = DateTimeField(null=True, blank=True)
    retrieval_end = DateTimeField(null=True, blank=True)
    retrieval_error = TextField(blank=True)

    class Meta:
        ordering = ('-created', 'dataset')
        unique_together = ('name', 'dataset')

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
            download = Download(self.url)
            download.save(self.dataset.file_directory, self.filename)
            self.retrieval_end = now()
        except Exception as error:
            self.retrieval_error = str(error)
            raise

        # We save the original filename instead of the new django one
        self.result = self.dataset.results.create(
            name=self.filename,
            fullpath=download.filepath,
            size=download.size,
        )
        self.save()

