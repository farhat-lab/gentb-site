import os
from hashlib import md5
import json

from model_utils.models import TimeStampedModel

from django.db import models
from django.core.files.storage import FileSystemStorage
from django.conf import settings
from django.utils.text import slugify
from django.core import serializers
from django.core.urlresolvers import reverse

from apps.tb_users.models import TBUser
from apps.utils.site_url_util import get_site_url

tb_file_system_storage = FileSystemStorage(location=settings.TB_SHARED_DATAFILE_DIRECTORY)

#def generate_new_filename(instance, filename):
#    #f, ext = os.path.splitext(filename)
#    instance.original_filename = basename(filename)
#    return join(instance.dataset.get_partial_path_for_datafile(), generate_storage_identifier())

DATASET_STATUS_UPLOADED_NOT_READY_ID = 1
DATASET_STATUS_UPLOADED_READY_ID = 2
DATASET_STATUS_PROCESSING_STARTED_ID = 3
DATASET_STATUS_PROCESSED_SUCCESS = 4
DATASET_STATUS_PROCESSED_FAILED = 4


class VCFDatasetStatus(models.Model):
    name = models.CharField(max_length=50)
    slug = models.SlugField(blank=True)
    sort_order = models.IntegerField()

    def __str__(self):
        return '%s - %s' % (self.sort_order, self.name)

    def save(self, *args, **kwargs):
        if not self.id:
            super(VCFDatasetStatus, self).save(*args, **kwargs)

        self.slug = slugify(self.name)

        super(VCFDatasetStatus, self).save(*args, **kwargs)

    class Meta:
        ordering = ('sort_order', '-name')
        verbose_name = 'VCF Dataset Status'
        verbose_name_plural = 'VCF Dataset Statuses'


class VCFDataset(TimeStampedModel):
    """
    Information from API call: https://api.github.com/repos/iqss/dataverse/milestones
    """
    user = models.ForeignKey(TBUser)

    title = models.CharField('Dataset title', max_length=255)

    description = models.TextField('Dataset description')

    status = models.ForeignKey(VCFDatasetStatus)

    file1 = models.FileField('File 1', upload_to='shared-files/%Y/%m'\
                            , storage=tb_file_system_storage)

    file2 = models.FileField('File 2 (optional)',
                             upload_to='shared-files/%Y/%m',
                             storage=tb_file_system_storage,
                             null=True,
                             blank=True)

    has_prediction = models.BooleanField('Has prediction results?',default=False, help_text='auto-filled on save')

    md5 = models.CharField(max_length=40, blank=True, db_index=True, help_text='auto-filled on save')


    def __str__(self):
        return self.title


    def run_script_link(self):
        if not self.id:
            return 'n/a'

        url = reverse('view_run_dataset_script', kwargs=dict(dataset_md5=self.md5))
        return '<a href="%s" target="_blank" style="display:block; border:1px solid #333; padding:10px; width:70px;">Run Script!</a>' % (url)
    run_script_link.allow_tags = True



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

    def save(self, *args, **kwargs):
        if not self.id:
            super(VCFDataset, self).save(*args, **kwargs)

        self.md5 = md5('%s%s' % (self.id, self.title)).hexdigest()

        super(VCFDataset, self).save(*args, **kwargs)

    def filename1(self):
        return os.path.basename(self.file1.name)

    def filename2(self):
        return os.path.basename(self.file2.name)


    def get_full_json(self):
        return serializers.serialize('json', VCFDataset.objects.filter(id=self.id))

    def get_script_args_json(self, run_md5, as_list=False):

        site_url = get_site_url()

        url_to_dataset = reverse('admin:predict_vcfdataset_change', args=(self.id,))
        admin_url = '{0}{1}'.format(site_url, url_to_dataset)
        callback_url = '{0}{1}'.format(site_url, reverse('view_dataset_run_notification', kwargs={}))

        d = dict(file1_path=self.file1.path,
                 dataset_id=self.id,
                 callback_url=callback_url,
                 user_email=self.user.user.email,
                 admin_url=admin_url,
                 run_md5=run_md5
                 )

        if self.file2:
            d['file2_path'] = self.file2.path

        if as_list:
            return [ json.dumps(d)]
            #return [ '\'%s\'' % json.dumps(d)]
        return json.dumps(d)

    def set_status_uploaded_ready(self, save_status=True):
        self.status = VCFDatasetStatus.objects.get(pk=DATASET_STATUS_UPLOADED_READY_ID)
        if save_status:
            self.save()

    def set_status_processing_started(self, save_status=True):
        self.status = VCFDatasetStatus.objects.get(pk=DATASET_STATUS_PROCESSING_STARTED_ID)
        if save_status:
            self.save()


    def set_status_process_failed(self, save_status=True):
        self.status = VCFDatasetStatus.objects.get(pk=DATASET_STATUS_PROCESSED_FAILED)
        self.has_prediction = False
        if save_status:
            self.save()

    def set_status_process_success(self, save_status=True):
        self.status = VCFDatasetStatus.objects.get(pk=DATASET_STATUS_PROCESSED_SUCCESS)
        self.has_prediction = True
        if save_status:
            self.save()


    class Meta:
        ordering = ('-created', 'title')
        verbose_name = 'VCF Dataset'
        verbose_name_plural = 'VCF Datasets'


class VCFDatasetNote(TimeStampedModel):

    dataset = models.ForeignKey(VCFDataset)

    title = models.CharField(max_length=255)
    note = models.TextField()

    class Meta:
        ordering = ('-modified', '-created')


class ScriptToRun(TimeStampedModel):

    name = models.CharField(max_length=100)
    is_chosen_script = models.BooleanField(default=True)
    script = models.TextField('Command line script run by webserver.  Arguments will be passed in JSON format.',
                              help_text='Example of JSON argument: \'{"admin_url": "http://127.0.0.1:8000/tb-admin/predict/vcfdataset/3/", "callback_url": "some_url to receive results", "dataset_id": 3, "user_email": "user_who_uploaded_file@place.edu", "file1_path": ".../tb_uploaded_files/shared-files/2015/08/Predict_-_genTB_BnVjFcO.png"}\'')

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

    dataset = models.ForeignKey(VCFDataset)

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

