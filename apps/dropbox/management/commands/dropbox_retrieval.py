"""
Examine existing PredictDataset objects.

- Retrieve files from confirmed dropbox links

"""

import logging
LOGGER = logging.getLogger('apps.dropbox')

from django.core.management.base import BaseCommand, CommandError
from apps.dropbox.models import DropboxFile

class Command(BaseCommand):
    help = """Regular run of new dropbox links:

      manage.py dropbox_retrieval

    """
    def handle(self, **options):
        print "Starting downloads..."
        for d_file in DropboxFile.objects.filter(result=''):
            print " * Downloading: %s" % d_file.url
            d_file.download_now()
            if d_file.retrieval_error:
                print "  XXX Failed: %s" % d_file.retrieval_error

