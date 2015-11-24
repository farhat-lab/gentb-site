import collections
from hashlib import md5
from django.utils import timezone

from model_utils.models import TimeStampedModel

from django.db import models
from django.core.urlresolvers import reverse

from jsonfield import JSONField # https://github.com/bradjasper/django-jsonfield

from apps.utils.site_url_util import get_site_url
from apps.dropbox_helper.forms import DropboxRetrievalParamsForm
from apps.predict.models import PredictDataset
from apps.utils.file_patterns import FilePatternHelper,\
                        FASTQ_PAIR_END_EXTENSION_TYPES

import logging
logger = logging.getLogger('apps.dropbox_helper.models')


class DropboxRetrievalLog(TimeStampedModel):

    dataset = models.OneToOneField(PredictDataset)

    # retrieved from dropbox
    file_metadata = JSONField(load_kwargs={'object_pairs_hook': collections.OrderedDict}, blank=True)
    file_metadata_err_msg = models.TextField(blank=True)

    # selected from metadata based on file endings
    selected_files = JSONField(load_kwargs={'object_pairs_hook': collections.OrderedDict}, blank=True)

    # system attempts to download files
    retrieval_start = models.DateTimeField(null=True, blank=True)
    retrieval_end = models.DateTimeField(null=True, blank=True)
    retrieval_error = models.TextField(blank=True)

    # Blank unless files are FastQ pair-end extensions
    fastq_pair_end_extension = models.CharField(max_length=20,\
                                    blank=True,\
                                    choices=FASTQ_PAIR_END_EXTENSION_TYPES,\
                                    help_text='For FastQ pair-end extensions. Either "_R" or "."')

    # success
    files_retrieved = models.BooleanField(default=False)

    # md5
    md5 = models.CharField(max_length=40, blank=True, db_index=True, help_text='auto-filled on save')

    def __str__(self):
        return '{0}'.format(self.dataset)

    def save(self, *args, **kwargs):
        if not self.id:
            super(DropboxRetrievalLog, self).save(*args, **kwargs)

        self.md5 = md5('%s%s' % (self.id, self.created)).hexdigest()

        self.fastq_pair_end_extension = self.get_fastq_extension_type()

        super(DropboxRetrievalLog, self).save(*args, **kwargs)

    def get_fastq_extension_type(self):
        """
        For pair-ended FastQ files, figure out the extension
        """
        return FilePatternHelper.get_fastq_extension_type(self.selected_files)

    class Meta:
        ordering = ('-created', 'dataset')
        #verbose_name = 'Dropbox Data Source'
        #verbose_name_plural = '{0}s'.format(verbose_name)

    def set_retrieval_start_time(self, with_reset=True):
        """
        Called before kicking off file retrieval.
        Note: This does not save the object
        """
        if with_reset:
            self.retrieval_end = None
            self.retrieval_error = ''
            self.files_retrieved = False

        self.retrieval_start = timezone.now()

    def set_retrieval_end_time(self, files_retrieved=True):
        """
        Called to mark the close of a file retrieval.
        Note: This does not save the object
        """
        if files_retrieved:
            self.files_retrieved = True

        self.retrieval_end = timezone.now()

    def get_dropbox_retrieval_script_params(self):
        assert self.id is not None, "This function cannot be called for an unsaved object.  (The id field is required)"

        callback_url = '{0}{1}'.format(get_site_url(),
                                reverse('record_file_retrieval_results', args=()))

        params = dict(dropbox_url=self.dataset.dropbox_url,
                destination_directory=self.dataset.file_directory,
                callback_url=callback_url,
                callback_md5=self.md5)

        print 'params', params

        f = DropboxRetrievalParamsForm(params)
        if f.is_valid():
            return f.cleaned_data


        logger.severe("Failed to create dropbox retrieve params for dataset id {0}. Errors from form: {1}".format(self.id, f.errors.items()))
        # Log a major error
        raise KeyError(f.errors)

"""
python manage.py shell
from apps.dropbox_helper.models import *
dlog = DropboxRetrievalLog.objects.first()
dlog.get_fastq_extension_type()

"""
