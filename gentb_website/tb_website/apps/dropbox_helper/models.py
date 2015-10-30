import collections
from hashlib import md5
from datetime import datetime
from model_utils.models import TimeStampedModel
from django.db import models
from django.core.urlresolvers import reverse
from jsonfield import JSONField # https://github.com/bradjasper/django-jsonfield
from apps.dropbox_helper.forms import DropboxRetrievalParamsForm
from apps.predict.models import PredictDataset

class DropboxRetrievalLog(TimeStampedModel):

    dataset = models.ForeignKey(PredictDataset, unique=True)

    # retrieved from dropbox
    file_metadata = JSONField(load_kwargs={'object_pairs_hook': collections.OrderedDict}, blank=True)
    file_metadata_err_msg = models.TextField(blank=True)

    # selected from metadata based on file endings
    selected_files = JSONField(load_kwargs={'object_pairs_hook': collections.OrderedDict}, blank=True)

    # system attempts to download files
    retrieval_start = models.DateTimeField(null=True, blank=True)
    retrieval_end = models.DateTimeField(null=True, blank=True)
    retrieval_error = models.TextField(blank=True)

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

        super(DropboxRetrievalLog, self).save(*args, **kwargs)

    class Meta:
        ordering = ('-created', 'dataset')
        #verbose_name = 'Dropbox Data Source'
        #verbose_name_plural = '{0}s'.format(verbose_name)

    def set_retrieval_start_time(self):
        self.retrieval_start = datetime.now()

    def set_retrieval_end_time(self):
        self.retrieval_end = datetime.now()

    def get_dropbox_retrieval_script_params(self):
        assert self.id is not None, "This function cannot be called for an unsaved object.  (The id field is required)"

        params = dict(dropbox_url=self.dataset.dropbox_url,
                destination_directory=self.dataset.file_directory,
                callback_url=reverse('record_file_retrieval_results', args=()),
                callback_md5=self.md5)

        f = DropboxRetrievalParamsForm(params)
        if f.is_valid():
            return f.cleaned_data

        # Log a major error
        return {}
