"""Download from urls any uploads from outside sources"""
import logging

from django.utils.timezone import now
from django.core.management.base import BaseCommand
from apps.uploads.models import DropboxUploadFile, ManualUploadFile, ResumableUploadFile

LOGGER = logging.getLogger('apps.uploads')

class Command(BaseCommand):
    """Run command for uploads"""
    help = __doc__ + """:

      manage.py process_uploads

    """
    @staticmethod
    def handle(**_):
        """Handle script call"""
        for cls in (DropboxUploadFile, ManualUploadFile):
            for d_file in cls.objects.filter(retrieval_start__isnull=True):
                print(" + Downloading: %s" % d_file.url)
                d_file.download_now()
                if d_file.retrieval_error:
                    print(" ! Error downloading")
        count = ResumableUploadFile.objects.filter(retrieval_start__isnull=True).update(
            retrieval_error='Retry resumable upload, can not happen in server.',
            retrieval_start=now(),
        )
        if count:
            print(" * Resumable uploads marked as impossible: {}".format(count))
