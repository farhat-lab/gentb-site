
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
    filename = CharField(max_length=255, null=True)
    size = PositiveIntegerField(null=True)
    icon = URLField(null=True)
    created = DateTimeField(auto_now_add=True)

    # system attempts to download files
    retrieval_start = DateTimeField(null=True, blank=True)
    retrieval_end = DateTimeField(null=True, blank=True)
    retrieval_error = TextField(blank=True)

    class Meta:
        ordering = ('-created', 'dataset')
        unique_together = ('filename', 'dataset')

    def __str__(self):
        return '{0} ({1})'.format(self.dataset, self.name)

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
            if not download.filepath:
                raise KeyError("No filename provided for download.")
        except Exception as error:
            self.retrieval_error = str(error)
            return

        # We save the original filename instead of the new django one
        self.result = self.dataset.results.create(
            name=self.filename,
            fullpath=download.filepath,
            size=download.size,
        )
        self.save()

