
import logging
LOGGER = logging.getLogger('apps.uploads')

from django.core.management.base import BaseCommand, CommandError
from apps.uploads.models import DropboxUploadFile, ManualUploadFile

class Command(BaseCommand):
    help = """Regular run of new dropbox links:

      manage.py process_uploads

    """
    def handle(self, **options):
        for cls in [DropboxUploadFile, ManualUploadFile]:
            for d_file in cls.objects.filter(retrieval_start__isnull=True):
                print " + Downloading: %s" % d_file.url
                d_file.download_now()
                if d_file.retrieval_error:
                    print " ! Error downloading"

