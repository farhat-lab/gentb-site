import os
import json

from hashlib import md5
from os.path import join, isfile, isdir

from collections import defaultdict
from model_utils.models import TimeStampedModel

from django.db import models
from django.conf import settings
from django.utils.text import slugify
from django.core import serializers
from django.core.urlresolvers import reverse

from apps.utils.site_url_util import get_site_url
from apps.utils.file_patterns import *

from apps.utils.result_file_info import RESULT_FILE_NAME_DICT,\
            EXPECTED_FILE_DESCRIPTIONS, RESULT_OUTPUT_DIRECTORY_NAME

from .script_runner import run_script

import logging
LOGGER = logging.getLogger('apps.predict.runner')

VCF_ANALYSIS_SCRIPT = 'analyseVCF.pl'
FASTQ_ANALYSIS_SCRIPT = 'analyseNGS.pl'
SCRIPT_DIR = join(settings.SITE_ROOT, 'apps', 'predict', 'predict_pipeline')


class PredictDatasetStatus(models.Model):
    name = models.CharField(max_length=50)
    human_name = models.CharField(max_length=100)
    is_error = models.BooleanField()
    slug = models.SlugField(blank=True)
    sort_order = models.IntegerField()

    def __str__(self):
        return '%s - %s' % (self.sort_order, self.name)

    def save(self, *args, **kwargs):
        if not self.pk:
            super(PredictDatasetStatus, self).save(*args, **kwargs)

        self.slug = slugify(self.name)

        super(PredictDatasetStatus, self).save(*args, **kwargs)

    class Meta:
        ordering = ('sort_order', '-name')
        #verbose_name = 'VCF Dataset Status'
        verbose_name_plural = 'Predict Dataset Statuses'


class PredictDataset(TimeStampedModel):
    """An uploaded predict dataset"""

    STATUS_DELETED = 0
    STATUS_NOT_READY = 1
    STATUS_CONFIRMED = 2
    STATUS_FILE_RETRIEVAL_STARTED = 3
    STATUS_FILE_RETRIEVAL_ERROR = 4
    STATUS_FILE_RETRIEVAL_COMPLETE = 5
    STATUS_PROCESSING_STARTED = 6
    STATUS_PROCESSED_SUCCESS = 7
    STATUS_PROCESSED_FAILED = 8

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='datasets')
    md5 = models.CharField(max_length=40, blank=True, db_index=True,
            help_text='auto-filled on save')
    title = models.CharField('Dataset title', max_length=255)
    file_type = models.CharField(choices=FILE_TYPES, max_length=25)

    fastq_type = models.CharField(max_length=50,\
        choices=FASTQ_FILE_TYPES, blank=True,\
        help_text='Only used for FastQ files')

    dropbox_url = models.URLField("Dropbox link", help_text='https://www.dropbox.com/help/274')
    description = models.TextField('Dataset description')
    status = models.ForeignKey(PredictDatasetStatus)
    file_directory = models.CharField(max_length=255, blank=True)
    has_prediction = models.BooleanField(default=False)

    def __str__(self):
        return self.title

    def get_absolute_url(self):
        return reverse('view_single_dataset', kwargs=dict(slug=self.md5))

    @property
    def files(self):
        return os.listdir(self.file_directory)

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

    def get_script_command(self):
        """
        Using dataset information, to decide whether to run:
            (1) script for a VCF file
            (2) script for FastQ files
        """
        # Formate either a VCF or FastQ pipeline command
        if self.is_vcf_file():     # (2a) command for a VCF file
            command_to_run = self.get_vcf_script_command()

        elif self.is_fastq_file(): # (2b) command for FastQ files
            command_to_run = self.get_fastq_script_command()

        else:
            err_title = 'Not VCF or FastQ file'
            err_note = 'Could not determine the file type.\
             Database contained: "%s"' % (self.file_type)
            #err_msg_obj =
            self.record_error(err_title, err_note)
            return None

        return command_to_run

    def get_pipeline_command(self):
        """
        Return the full pipeline command to run the appropriate analyze script
        for this dataset's files
        """
        try:
            return (True, self.get_script_command())
        except Exception as err:
            return (False, err.args[0])

    def get_fastq_script_command(self):
        """
        Run the command for FastQ files.
        e.g. perl analyseNGS.pl (directory name)
        """
        if self.err_found:
            return None

        # (1) Make sure the 'analyseNGS.pl' command is in the
        #   specified 'Pipeline Scripts Directory'
        #
        script_cmd = join(SCRIPT_DIR, FASTQ_ANALYSIS_SCRIPT)
        if not isfile(script_cmd):
            err_title = 'The "%s" file was not found' % (FASTQ_ANALYSIS_SCRIPT)
            err_note = """The "%s" file was not found at this location: \n%s\n
            To change the path to the script, please go into the admin control\
             panel and modify the 'Pipeline Scripts Directory'.\n\
            Talk to your administrator for details.""" % (FASTQ_ANALYSIS_SCRIPT, script_cmd)
            #err_msg_obj = self.record_error(err_title, err_note)
            self.record_error(err_title, err_note)
            return None

        # Format the full command with target containing
        #   input files
        #
        if self.is_fastq_single_ended():
            command_str = 'perl {0} 0 . {1}'.format(script_cmd,\
                self.file_directory)
        elif self.is_fastq_pair_ended():
            pair_extension = self.get_fastq_pair_end_extension()
            if pair_extension is None:
                err_title = 'FastQ could not find pair-ended extension type'
                err_note = 'Could not determine pair-ended extension type.\
                 Database contained: "%s"' % (self.fastq_type)
                #err_msg_obj = self.record_error(err_title, err_note)
                self.record_error(err_title, err_note)
                return None

            command_str = 'perl {0} 1 {1} {2}'.format(script_cmd,\
                             pair_extension, self.file_directory)
        else:
            err_title = 'FastQ: single-ended or pair-ended?'
            err_note = 'Could not determine single-ended or pair-ended FastQ type.\
             Database contained: "%s"' % (self.fastq_type)
            # err_msg_obj = self.record_error(err_title, err_note)
            self.record_error(err_title, err_note)
            return None

        return command_str

    def get_vcf_script_command(self):
        """
        Run the command for a VCF file.
        e.g. perl analyseVCF.pl (directory name)
        """
        # (1) Make sure the 'analyseVCF.pl' command is in the
        #   specified 'Pipeline Scripts Directory'
        #
        script_cmd = join(SCRIPT_DIR, VCF_ANALYSIS_SCRIPT)
        if not isfile(script_cmd):
            err_title = 'The "%s" file was not found' % (VCF_ANALYSIS_SCRIPT)
            err_note = """The "%s" file was not found at this location: \n%s\n
            To change the path to the script, please go into the admin control\
             panel and modify the 'Pipeline Scripts Directory'.\n\
            Talk to your administrator for details.""" % (VCF_ANALYSIS_SCRIPT, script_cmd)
            #err_msg_obj = self.record_error(err_title, err_note)
            self.record_error(err_title, err_note)
            return None

        # Format the full command with target containing
        #   input files
        #
        command_str = 'perl {0} {1}'.format(script_cmd,\
                self.file_directory)

        return command_str

    def run_command(self):
        """
        - (1) Create a DatasetScriptRun object
        - (2) Make a custom callback script add the dataset file directory
        - (3) Update the Dataset status
        - (4) Run the script
        """
        ready, command_to_run = self.get_pipeline_command()
        if not ready:
            return (False, command_to_run)

        # (1) Create a run object and save the command being run
        #
        dsr = DatasetScriptRun(dataset=self)
        dsr.notes = command_to_run
        dsr.save()

        # ---------------------------------------------------------
        # (2) Place a callback script with the Dataset directory
        #   - contains callback url + md5
        #   - called after the analyse scripts are Run
        # ---------------------------------------------------------

        # Get callback args
        template_dict = dict(callback_info_dict=dsr.callback_info())
        template_dict.update(RESULT_FILE_NAME_DICT)
        template_dict['RESULT_OUTPUT_DIRECTORY_NAME'] = RESULT_OUTPUT_DIRECTORY_NAME
        template_dict['EXPECTED_FILE_DESCRIPTIONS'] = EXPECTED_FILE_DESCRIPTIONS
        

        print('-' * 40)
        print(template_dict)
        print('-' * 40)

        # Place script in dataset file_directory
        script_fullname = join(self.file_directory, 'feedback.json')
        fhandler = open(script_fullname, 'w')
        fhandler.write(json.dumps(template_dict))
        fhandler.close()

        # (3) Update the dataset status to 'in process'
        #
        self.set_status_processing_started()

        # (4) Run the script -- not waiting for output
        #
        full_args = command_to_run.split()
        print ('full_args', full_args)
        run_script(full_args)

    def record_error(self, msg_title, msg):
        """
        Record an error:
            - within this non-persistnt object
            - in the log files
            - in the database in a PredictDatasetNote object
        """
        # Write to the logs
        #
        LOGGER.error(msg_title)
        LOGGER.error(msg)

        # Write to the database
        self.notes.create(title=msg_title, note=msg)

    def get_fastq_pair_end_extension(self):
        try:
            dlog = self.dropboxretrievallog
        except:
            return None

        return dlog.fastq_pair_end_extension

    def get_file_patterns(self):
        return FilePatternHelper.get_file_patterns_for_dropbox(self.file_type)

    def get_heatmap(self):
        data = None
        maf = os.path.join(self.file_directory, 'output', 'matrix.json')
        if not os.path.isfile(maf):
            return None

        ret = defaultdict(list)
        with open(maf, 'r') as fhl:
            data = json.loads(fhl.read())
        for row in data[0]:
            ret['data'].append(row[2])
            if row[0] not in ret['rows']:
                ret['rows'].append(row[0])
            if row[1] not in ret['cols']:
                ret['cols'].append(row[1])
        ret['dim'] = [len(ret['rows']), len(ret['cols'])]
        return ret

    def user_name(self):
        if self.user:
            return self.user.username
        return 'n/a'

    def user_affiliation(self):
        if self.user:
            return self.user.affiliation
        return 'n/a'

    def user_email(self):
        if self.user:
            return self.user.email
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
        self.set_status(self.STATUS_NOT_READY, save_status)

    def set_status_confirmed(self, save_status=True):
        self.set_status(self.STATUS_CONFIRMED, save_status)

    # File Retrieval statuses
    #
    def set_status_file_retrieval_started(self, save_status=True):
        self.set_status(self.STATUS_FILE_RETRIEVAL_STARTED, save_status)

    def set_status_file_retrieval_error(self, save_status=True):
        self.set_status(self.STATUS_FILE_RETRIEVAL_ERROR, save_status)

    def set_status_file_retrieval_complete(self, save_status=True):
        self.set_status(self.STATUS_FILE_RETRIEVAL_COMPLETE, save_status)

    # Pipeline processing statuses
    #
    def set_status_processing_started(self, save_status=True):
        self.set_status(self.STATUS_PROCESSING_STARTED, save_status)

    def set_status_processing_success(self, save_status=True):
        self.set_status(self.STATUS_PROCESSED_SUCCESS, save_status)

    def set_status_processing_failed(self, save_status=True):
        self.set_status(self.STATUS_PROCESSED_FAILED, save_status)


    class Meta:
        ordering = ('-created', 'title')


class PredictDatasetNote(TimeStampedModel):
    """Notes of background processes"""
    dataset = models.ForeignKey(PredictDataset, related_name='notes')
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
    dataset = models.ForeignKey(PredictDataset, related_name='runs')
    md5 = models.CharField(max_length=40, blank=True, db_index=True)

    notes = models.TextField(blank=True)
    result_received = models.BooleanField(default=False)
    result_success = models.BooleanField(default=False)
    result_data = models.TextField(blank=True)

    def __str__(self):
        return '%s' % self.dataset

    def callback_info(self):
        """
        Creates a dictionary of data to send which links back to us here.
        """
        pk = self.dataset.pk
        admin_url = reverse('admin:predict_predictdataset_change', args=[pk])

        slug = {'slug': self.md5}
        callback_url = get_site_url(for_internal_callback=True) + \
                reverse('view_dataset_run_notification', kwargs=slug)

        return dict(
             dataset_id=pk,
             file_directory=self.dataset.file_directory,
             callback_url=callback_url,
             user_email=self.dataset.user.email,
             admin_url=get_site_url() + admin_url,
             run_md5=self.md5
         )

    def save(self, *args, **kwargs):
        if not self.md5:
            self.md5 = md5('%s%s' % (self.dataset.pk, self.created)).hexdigest()

        if self.result_received:
            if self.result_success:
                self.dataset.set_status_processing_success()
            else:
                self.dataset.set_status_processing_failed()

        super(DatasetScriptRun, self).save(*args, **kwargs)

    class Meta:
        ordering = ('-modified', '-created')
