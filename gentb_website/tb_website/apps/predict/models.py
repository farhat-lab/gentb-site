import collections
from hashlib import md5
import json
from os.path import basename, join, isdir
import os

from datetime import datetime

from model_utils.models import TimeStampedModel

from django.db import models
from django.core.files.storage import FileSystemStorage
from django.conf import settings
from django.utils.text import slugify
from django.core import serializers
from django.core.urlresolvers import reverse
from jsonfield import JSONField # https://github.com/bradjasper/django-jsonfield

from apps.tb_users.models import TBUser
from apps.utils.site_url_util import get_site_url
from apps.utils.file_patterns import *

#tb_file_system_storage = FileSystemStorage(location=settings.TB_SHARED_DATAFILE_DIRECTORY)

#def generate_new_filename(instance, filename):
#    #f, ext = os.path.splitext(filename)
#    instance.original_filename = basename(filename)
#    return join(instance.dataset.get_partial_path_for_datafile(), generate_storage_identifier())

DATASET_STATUS_NOT_READY_ID = 1
DATASET_STATUS_CONFIRMED_ID = 2

DATASET_STATUS_FILE_RETRIEVAL_STARTED = 3
DATASET_STATUS_FILE_RETRIEVAL_ERROR = 4
DATASET_STATUS_FILE_RETRIEVAL_COMPLETE = 5

DATASET_STATUS_PROCESSING_STARTED_ID = 6
DATASET_STATUS_PROCESSED_SUCCESS = 7
DATASET_STATUS_PROCESSED_FAILED = 8



class PredictDatasetStatus(models.Model):
    name = models.CharField(max_length=50)
    human_name = models.CharField(max_length=100)
    is_error = models.BooleanField()
    slug = models.SlugField(blank=True)
    sort_order = models.IntegerField()

    def __str__(self):
        return '%s - %s' % (self.sort_order, self.name)

    def save(self, *args, **kwargs):
        if not self.id:
            super(PredictDatasetStatus, self).save(*args, **kwargs)

        self.slug = slugify(self.name)

        super(PredictDatasetStatus, self).save(*args, **kwargs)

    class Meta:
        ordering = ('sort_order', '-name')
        #verbose_name = 'VCF Dataset Status'
        verbose_name_plural = 'Predict Dataset Statuses'


class PredictDataset(TimeStampedModel):
    """

    """
    user = models.ForeignKey(TBUser)

    title = models.CharField('Dataset title', max_length=255)

    file_type = models.CharField(choices=FILE_TYPES, max_length=25)

    fastq_type = models.CharField(max_length=50,\
        choices=FASTQ_FILE_TYPES,\
        blank=True,\
        help_text='Only used for FastQ files')

    dropbox_url = models.URLField("Dropbox link", help_text='https://www.dropbox.com/help/274')

    description = models.TextField('Dataset description')

    status = models.ForeignKey(PredictDatasetStatus)

    file_directory = models.CharField(max_length=255, blank=True)

    has_prediction = models.BooleanField(default=False)

    md5 = models.CharField(max_length=40, blank=True, db_index=True, help_text='auto-filled on save')

    def __str__(self):
        return self.title


    def is_vcf_file(self):
        return FilePatternHelper.is_vcf_file(self.file_type)

    def is_fastq_file(self):
        return FilePatternHelper.is_fastq_file(self.file_type)

    def is_fastq_single_ended(self):
        if not self.is_fastq_file():
            return False
        return FilePatternHelper.is_fastq_single_ended(self.fastq_type)

    def is_fastq_pair_ended(self):
        if not self.is_fastq_file():
            return False
        return FilePatternHelper.is_fastq_pair_ended(self.fastq_type)

    def get_fastq_pair_end_extension(self):
        try:
            dlog = self.dropboxretrievallog
        except:
            return None

        return dlog.fastq_pair_end_extension

    def get_file_patterns(self):
        return FilePatternHelper.get_file_patterns_for_dropbox(self.file_type)


    def user_name(self):
        if self.user:
            return self.user.user
        return 'n/a'

    def user_affiliation(self):
        if self.user:
            return self.user.affiliation
        return 'n/a'

    def user_email(self):
        if self.user:
            return self.user.user.email
        return 'n/a'
    #'user_name', 'user_email'

    def create_dataset_directory_name(self):
        """
        Create a directory based on the id of this object.

        e.g. id = 5
        dirname = 'tbdata_00000005'
        Attempt to create the directory under settings.TB_SHARED_DATAFILE_DIRECTORY
        """
        assert self.id is not None, "The object must be saved (and have an 'id') before this method is called."

        # zero pad the object id
        #
        job_num = str(self.id).zfill(8)

        # directory name is (tb_file_system_storage + "tb_data_" + job_num)
        #
        dirname = join(settings.TB_SHARED_DATAFILE_DIRECTORY, 'tbdata_{0}'.format(job_num))

        # create the new directory (if it doesn't exist)
        if not isdir(dirname):
            os.makedirs(dirname)

        return dirname

        #datetime.now().strftime('%Y-%m-%d_%H-%M-%S')



    def save(self, *args, **kwargs):
        if not self.id:
            super(PredictDataset, self).save(*args, **kwargs)

        # Set the md5
        self.md5 = md5('%s%s' % (self.id, self.title)).hexdigest()

        # -----------------------------
        # Initialize the file directory
        # for this dataset
        # -----------------------------
        if not self.file_directory:
            self.file_directory = self.create_dataset_directory_name()

        super(PredictDataset, self).save(*args, **kwargs)


    def get_full_json(self):
        """
        Need to serialize the PredictDataset as well as related PredictDatasetFile objects
        """
        return serializers.serialize('json', PredictDataset.objects.filter(id=self.id))


    def get_script_args_json(self, run_md5, **kwargs):

        as_list = kwargs.get('as_list', False)
        as_dict = kwargs.get('as_dict', False)

        site_url = get_site_url()

        url_to_dataset = reverse('admin:predict_predictdataset_change', args=(self.id,))
        admin_url = '{0}{1}'.format(site_url, url_to_dataset)
        callback_url = '{0}{1}'.format(site_url, reverse('view_dataset_run_notification', kwargs={}))

        d = dict(file_directory=self.file_directory,
                 dataset_id=self.id,
                 callback_url=callback_url,
                 user_email=self.user.user.email,
                 admin_url=admin_url,
                 run_md5=run_md5
                 )

        if as_dict:
            return d

        if as_list:
            return [ json.dumps(d)]
            #return [ '\'%s\'' % json.dumps(d)]
        return json.dumps(d)

    def set_status(self, status_type, save_status=True):
        try:
            new_status = PredictDatasetStatus.objects.get(pk=status_type)
        except PredictDatasetStatus.DoesNotExist:
            return

        self.status = new_status
        if save_status:
            self.save()

    # Initial information statuses
    #
    def set_status_not_ready(self, save_status=True):
        self.set_status(DATASET_STATUS_NOT_READY_ID, save_status)

    def set_status_confirmed(self, save_status=True):
        self.set_status(DATASET_STATUS_CONFIRMED_ID, save_status)

    # File Retrieval statuses
    #
    def set_status_file_retrieval_started(self, save_status=True):
        self.set_status(DATASET_STATUS_FILE_RETRIEVAL_STARTED, save_status)

    def set_status_file_retrieval_error(self, save_status=True):
        self.set_status(DATASET_STATUS_FILE_RETRIEVAL_ERROR, save_status)

    def set_status_file_retrieval_complete(self, save_status=True):
        self.set_status(DATASET_STATUS_FILE_RETRIEVAL_COMPLETE, save_status)

    # Pipeline processing statuses
    #
    def set_status_processing_started(self, save_status=True):
        self.set_status(DATASET_STATUS_PROCESSING_STARTED_ID, save_status)

    def set_status_processing_success(self, save_status=True):
        self.set_status(DATASET_STATUS_PROCESSED_SUCCESS, save_status)

    def set_status_processing_failed(self, save_status=True):
        self.set_status(DATASET_STATUS_PROCESSED_FAILED, save_status)


    class Meta:
        ordering = ('-created', 'title')


class PredictDatasetNote(TimeStampedModel):

    dataset = models.ForeignKey(PredictDataset)

    title = models.CharField(max_length=255)
    note = models.TextField()

    def __str__(self):
        return self.title

    class Meta:
        ordering = ('-modified', '-created')

class PredictDatasetFile(TimeStampedModel):
    """
    Single file that has been retrieved for a particular dataset
    """
    dataset = models.ForeignKey(PredictDataset)
    name = models.CharField(max_length=255, help_text='Name of the file (w/o) the path')
    fullpath = models.TextField(help_text='Full path to the file')
    size = models.IntegerField(default=0, help_text='Size of the file in bytes')

    def __str__(self):
        return self.name

    class Meta:
        ordering = ('-created', 'dataset', 'name')

class PipelineScriptsDirectory(TimeStampedModel):
    """
    Give the directory containinng Perls scripts: analyseNGS.pl and analyseVCF.pl
    """
    name = models.CharField(max_length=100, help_text='helpful user name')
    script_directory = models.TextField(help_text='Full path to the directory \
    containing the analyseNGS.pl and analyseVCF.pl pipeline scripts')
    is_chosen_directory = models.BooleanField(default=True)

    class Meta:
        ordering = ('is_chosen_directory', '-modified')
        verbose_name = 'Pipeline Scripts Directory'
        verbose_name_plural = verbose_name

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):

        # Strip the directory name and, if needed, remove ending file separator
        self.script_directory = self.script_directory.strip()

        super(PipelineScriptsDirectory, self).save(*args, **kwargs)

class ScriptToRun(TimeStampedModel):

    name = models.CharField(max_length=100)
    is_chosen_script = models.BooleanField(default=True)
    script = models.TextField('Command line script run by webserver.  Arguments will be passed in JSON format.',
                              help_text='Example of JSON argument: \'{"admin_url": "http://127.0.0.1:8000/tb-admin/predict/PredictDataset/3/", "callback_url": "some_url to receive results", "dataset_id": 3, "user_email": "user_who_uploaded_file@place.edu", "file1_path": ".../tb_uploaded_files/shared-files/2015/08/Predict_-_genTB_BnVjFcO.png"}\'')

    script_args = models.TextField(blank=True, help_text='populated on save')

    def __str__(self):
        return self.name

    def get_script_args_as_list(self):
        return self.script.split()

    def save(self, *args, **kwargs):

        self.script_args = self.script.split()

        try:
            ot = ScriptToRun.objects.get(is_chosen_script=True)
            if not self==ot:
                ot.is_chosen_script = False
                ot.save()
        except ScriptToRun.DoesNotExist:
            pass

        super(ScriptToRun, self).save(*args, **kwargs)

class DatasetScriptRun(TimeStampedModel):

    dataset = models.ForeignKey(PredictDataset)

    notes = models.TextField(blank=True)

    result_received = models.BooleanField(default=False)

    result_success = models.BooleanField(default=False)

    result_data = models.TextField(blank=True)

    md5 = models.CharField(max_length=40, blank=True, db_index=True, help_text='auto-filled on save')

    def __str__(self):
        return '%s' % self.dataset

    def save(self, *args, **kwargs):
        if not self.id:
            super(DatasetScriptRun, self).save(*args, **kwargs)

        self.md5 = md5('%s%s' % (self.id, self.created)).hexdigest()

        super(DatasetScriptRun, self).save(*args, **kwargs)

    class Meta:
        ordering = ('-modified', '-created')
