import os
import json

from hashlib import md5
from os.path import join, isfile, isdir, basename

from collections import defaultdict
from model_utils.models import TimeStampedModel

from django.db import models
from django.conf import settings
from django.utils.text import slugify
from django.core import serializers
from django.core.urlresolvers import reverse
from django.utils.timezone import now

from django.utils.translation import ugettext_lazy as _
from apps.mutations.models import Drug
from apps.utils.result_file_info import RESULT_FILE_NAME_DICT,\
            EXPECTED_FILE_DESCRIPTIONS, RESULT_OUTPUT_DIRECTORY_NAME

from .script_runner import run_script
from .utils import *

import logging
LOGGER = logging.getLogger('apps.predict.pipeline')

class PredictDataset(TimeStampedModel):
    """An uploaded predict dataset"""
    # These names are also keys, DO NOT CHANGE
    STATUS_CHOICES = list(enumerate([
      _('Dataset Deleted'),
      _('Dataset Not Ready'),
      _('Dataset Confirmed'),
      _('File Retrieval Started'),
      _('File Retrieval Failed'),
      _('File Retrieval Success'),
      _('Processing Started'),
      _('Processing Success'),
      _('Processing Failed'),
    ]))
    STATUS = dict([(x, x) for x, st in STATUS_CHOICES])
    STATUS.update(dict([(str(st).upper().replace(' ', '_'), x) for x, st in STATUS_CHOICES]))
    STATUS_ERRORS = [0, 5, 9]

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='datasets')
    md5 = models.CharField(max_length=40, blank=True, db_index=True,
            help_text='auto-filled on save')
    title = models.CharField('Dataset Title', max_length=255)
    file_type = models.CharField(choices=FILE_TYPES, max_length=25)

    fastq_type = models.CharField(max_length=50,\
        choices=FASTQ_FILE_TYPES, blank=True,\
        help_text='Only used for FastQ files')

    description = models.TextField('Dataset description')
    file_directory = models.CharField(max_length=255, blank=True)
    has_prediction = models.BooleanField(default=False)

    status = models.PositiveIntegerField(default=1, choices=STATUS_CHOICES)
    is_error = property(lambda self: self.status in self.STATUS_ERRORS)

    def __str__(self):
        return self.title

    def get_absolute_url(self):
        return reverse('predict:view_single_dataset', kwargs=dict(slug=self.md5))

    @property
    def result_files(self):
        try:
            return os.listdir(join(self.file_directory, 'output'))
        except OSError:
            return []

    @property
    def media_url(self):
        return join(settings.MEDIA_URL, 'data', basename(self.file_directory))

    def check_for_prediction(self):
        if isfile(join(self.file_directory, 'output', 'matrix.json')):
            self.has_prediction = True
            self.save()

    def is_manual(self):
        return self.file_type == 'manual'

    def get_script_command(self):
        """
        Using dataset information, to decide whether to run:
            (1) script for a VCF file
            (2) script for FastQ files
            (3) script for Manual input
        """
        if self.file_type == FILE_TYPE_VCF:
            command_to_run = self.get_vcf_script_command()
        elif self.file_type == FILE_TYPE_FASTQ:
            command_to_run = self.get_fastq_script_command()
        elif self.file_type == FILE_TYPE_MANUAL:
            command_to_run = self.get_manual_script_command()

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

        return ' '.join(['perl', script_cmd,
            str(int(self.fastq_type == FASTQ_PAIR_ENDED)),
            FASTQ_PAIR_END[self.fastq_type],
            self.file_directory])

    def get_manual_script_command(self):
        """Manual script processing"""
        script_cmd = join(SCRIPT_DIR, MANUAL_ANALYSIS_SCRIPT)
        return 'perl {0} {1}'.format(script_cmd, self.file_directory)

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
        self.set_status('PROCESSING_STARTED')

        # (4) Run the script -- not waiting for output
        #
        full_args = command_to_run.split()
        print ('full_args', full_args)
        run_script(full_args)
        return (True, "DONE")

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

    def make_scatter(self, locusts, data):
        regions = defaultdict(list)
        for gene in data:
            if gene:
                region = gene.split("_")[-1].lower()
                regions[region].append(gene)

        for x, locust in enumerate(locusts):
            key = locust.lower()
            ret = {"x": x, "y": 0, "size": 5, "tip": ["No mutations"]}
            if key in regions:
                ret["y"] = len(regions[key])
                ret["size"] = 9
                ret["tip"] = regions[key]
            yield ret

    def get_heatmap(self):
        data = None
        maf = os.path.join(self.file_directory, 'output', 'matrix.json')
        if not os.path.isfile(maf):
            return None

        ret = defaultdict(list)
        drugs = list()
        with open(maf, 'r') as fhl:
            data = json.loads(fhl.read())
        for row in data[0]:
            ret['data'].append(row[2])
            ret['extra'].append([row[3], row[4]])
            if row[0] not in ret['rows']:
                ret['rows'].append(row[0])
            if row[1] not in ret['cols']:
                ret['cols'].append(row[1])

            try:
                drugs.append(Drug.objects.get(code__iexact=row[1]))
            except Drug.DoesNotExist:
                import sys
                sys.stderr.write("Can't find drug %s\n" % row[1])
                drugs.append(None)

        ret['dim'] = [len(ret['rows']), len(ret['cols'])]

        output = defaultdict(lambda: defaultdict(lambda: [None] * 2))
        for series, rows in enumerate(data[1:]):
            for row in rows:
                for col, datum in enumerate(zip(*rows[row])):
                    drug = drugs[col]
                    if drug is None:
                        output[row][col][series] = {}
                        continue
                    
                    locusts = list(drug.gene_locuses.values_list('name', flat=True))
                    output[row][col][series] = {
                        "cols": locusts,
                        "key": ["Important", "Other"][series],
                        "color": ["rgba(255, 0, 0, 0.8)", "rgba(0, 0, 255, 0.17)"][series],
                        "yAxis": "1",
                        "values": list(self.make_scatter(locusts, datum)),
                    }

        ret['scatter'] = {
          'data': output,
        }
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
        Need to serialize the PredictDataset
        """
        return serializers.serialize('json', PredictDataset.objects.filter(id=self.id))

    def set_status(self, status, save=True):
        """Set the status of this Dataset.
        
           status - Can be integer (0-9) or a tag such as 'FILE_RETRIEVAL_SUCCESS'
           save   - Boolean, if save should be called (default=True)
        """
        try:
            self.status = self.STATUS[status]
            if save:
                self.save()
        except KeyError:
            raise ValueError("Status %s not acceptable choice" % str(status))

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
    dataset = models.ForeignKey(PredictDataset, related_name='results')
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
    script = models.TextField('Command Line Script',
        help_text='Example of JSON argument: \'{"admin_url": "http://127.0.'
        '0.1:8000/tb-admin/predict/PredictDataset/3/", "callback_url": "som'
        'e_url to receive results", "dataset_id": 3, "user_email": "user_wh'
        'o_uploaded_file@place.edu", "file1_path": ".../tb_uploaded_files/s'
        'hared-files/2015/08/Predict_-_genTB_BnVjFcO.png"}\'')
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

    process_start = models.DateTimeField(auto_now_add=True, null=True)
    process_end = models.DateTimeField(null=True, blank=True)

    def __str__(self):
        return '%s' % self.dataset

    def callback_info(self):
        """
        Creates a dictionary of data to send which links back to us here.
        """
        pk = self.dataset.pk
        admin_url = reverse('admin:predict_predictdataset_change', args=[pk])

        slug = {'slug': self.md5}
        callback_url = get_site_url(internal=True) + \
                reverse('predict:callback', kwargs=slug)

        return dict(
             dataset_id=pk,
             file_directory=self.dataset.file_directory,
             callback_url=callback_url,
             user_email=self.dataset.user.email,
             admin_url=get_site_url() + admin_url,
             run_md5=self.md5
         )

    @property
    def process_time(self):
        if self.process_end and self.process_start:
            return self.process_end - self.process_start

    def save(self, *args, **kwargs):
        if not self.md5:
            self.md5 = md5('%s%s' % (self.dataset.pk, self.created)).hexdigest()

        if self.result_received:
            if not self.process_end:
                self.process_end = now()
            if self.result_success:
                self.dataset.check_for_prediction()
                self.dataset.set_status('PROCESSING_SUCCESS')
            else:
                self.dataset.set_status('PROCESSING_FAILED')

        super(DatasetScriptRun, self).save(*args, **kwargs)

    class Meta:
        ordering = ('-modified', '-created')
