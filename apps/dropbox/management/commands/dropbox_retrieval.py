"""
Examine existing PredictDataset objects.

- Retrieve files from confirmed dropbox links

"""

import logging
LOGGER = logging.getLogger('apps.dropbox')

from django.core.management.base import BaseCommand, CommandError

from apps.predict.models import PredictDataset
from apps.dropbox.models import DropboxFile

class RetrievalError(IOError):
    pass

class Command(BaseCommand):
    help = """Regular run of new dropbox links:

      manage.py dropbox_retrieval

    """
    def handle(self, **options):
        for dataset in PredictDataset.objects.filter(status=PredictDataset.STATUS_CONFIRMED):
            print " * Looking at dataset: %s" % str(dataset)
            dataset.set_status_file_retrieval_started()
            try:
                for d_file in DropboxFile.objects.filter(result__isnull=True):
                    print "  + Downloading: %s" % d_file.url
                    d_file.download_now()
                    if d_file.retrieval_error:
                        raise RetrievalError(d_file.retrieval_error)
            except RetrievalError as err:
                print "  XXX Failed: %s" % str(err)
                dataset.set_status_file_retrieval_error()
            else:
                dataset.set_status_file_retrieval_complete()

