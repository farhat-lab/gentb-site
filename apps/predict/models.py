#
# Copyright (C) 2017 Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Provide prediction app using the pipeline for building predictions and the
uploads app to download large data files from the users.
"""

import sys
import json
import logging

from hashlib import md5

import os
from os.path import join, isdir, basename

from collections import defaultdict
from datetime import timedelta

from model_utils.models import TimeStampedModel

from django.db.models import (
    Model, SET_NULL,
    ForeignKey, CharField, SlugField, TextField, BooleanField,
)
from django.conf import settings
from django.utils.text import slugify
from django.core import serializers
from django.core.urlresolvers import reverse
from django.utils.timezone import now
from django.utils.translation import ugettext_lazy as _

from apps.uploads.models import UploadFile
from apps.pipeline.models import Pipeline, PipelineRun, ProgramRun
from apps.mutations.models import Drug
from apps.mutations.utils import unpack_mutation_format

LOGGER = logging.getLogger('apps.predict')

# Basic status matrix, file_status + processing_status
STATUS_LABELS = [
    _('Dataset Confirmed'),       # STATUS_WAIT
    _('File Retrieval Started'),  # STATUS_START
    _('File Retrieval Failed'),   # STATUS_ERROR
    _('File Retrieval Success'),  # STATUS_DONE + STATUS_WAIT
    _('Processing Started'),      # STATUS_DONE + STATUS_START
    _('Processing Failed'),       # STATUS_DONE + STATUS_ERROR
    _('Processing Success'),      # STATUS_DONE + STATUS_DONE
    _('Prediction Ready'),        # STATUS_READY
    _('No Strains to process'),   # STATUS_NONE
    _('Processing Timed Out'),    # STATUS_TIMEOUT
    _('No Files Uploaded'),       # STATUS_NOFILES
]
(STATUS_WAIT, STATUS_START, STATUS_ERROR, STATUS_DONE) = range(4)
(STATUS_READY, STATUS_NONE) = (7, 8)
(STATUS_TIMEOUT, STATUS_NOFILES) = (-2, -1)
STATUS_LEVELS = (['default', 'primary', 'danger'] * 2) + \
    ['default', 'success', 'danger', 'danger', 'danger']

def get_timeout(timeout=14):
    """Returns the timedate when the prediction should time out"""
    return now() - timedelta(days=timeout)

class PredictDataset(TimeStampedModel):
    """An uploaded predict dataset"""
    BASE_DIR = settings.TB_SHARED_DATAFILE_DIRECTORY

    FILE_TYPE_VCF = 'vcf'
    FILE_TYPE_FASTQ = 'fastq'
    FILE_TYPE_FASTQ2 = 'fastq-pair'
    FILE_TYPE_MANUAL = 'manual'
    FILE_TYPES = [
        (FILE_TYPE_VCF, 'Variant Call Format (VCF)'),
        (FILE_TYPE_FASTQ, 'FastQ Single Ended Nucleotide Sequence'),
        (FILE_TYPE_FASTQ2, 'FastQ Pair Ended Nucleotide Sequences'),
        (FILE_TYPE_MANUAL, 'Mutations Manual Entry'),
    ]

    user = ForeignKey(settings.AUTH_USER_MODEL, related_name='datasets')
    md5 = CharField(max_length=40, blank=True, db_index=True,\
            help_text='auto-filled on save')
    title = CharField('Dataset Title', max_length=255)
    file_type = CharField(choices=FILE_TYPES, max_length=25)
    delete_sources = BooleanField(default=False,\
        help_text="If this is checked, we will delete all your input files"\
            " downloaded from dropbox after running the predict.")

    description = TextField('Dataset description')
    file_directory = CharField(max_length=255, blank=True)

    class Meta:
        ordering = ('-created', 'title')

    def __str__(self):
        return str(self.title)

    def get_absolute_url(self):
        """Return a link to thedataset view"""
        return reverse('predict:view_single_dataset', kwargs=dict(slug=self.md5))

    def get_status(self):
        """Returns the status as a string"""
        return STATUS_LABELS[self.status]

    def get_status_level(self):
        """Returns the btn/bootstrap color level for this status"""
        return STATUS_LEVELS[self.status]

    @property
    def status(self):
        """Returns the numeric level which this status has"""
        return min([strain.status for strain in self.strains.all()] + [STATUS_NONE])

    @property
    def directory_exists(self):
        """Returns true if the file_directory exists"""
        return os.path.isdir(self.file_directory)

    @property
    def has_prediction(self):
        """Returns true if any strain has a predict file"""
        return any([strain.has_prediction for strain in self.strains.all()])

    @property
    def has_lineages(self):
        """Return true if any strain has a lineage file"""
        return any([strain.has_lineage for strain in self.strains.all()])

    @property
    def has_output_files(self):
        """Return true if any strain has output files"""
        return any([list(strain.output_files) for strain in self.strains.all()])

    @property
    def time_taken(self):
        strain = self.strains.first()

        # check if these exist
        if not strain:
            return None

        if not strain.piperun:
            return 'Error: pipeline not started'

        prs = ProgramRun.objects.filter(piperun=strain.piperun)

        # one of the program runs had an error
        if prs.filter(is_error=True).count():
            return 'Error in pipeline'

        # time of the final program run completion or now if still running
        dt = now()
        if not prs.filter(completed__isnull=True).count():
            pr = prs.filter(completed__isnull=False).order_by('-completed').first()
            if pr:
                dt = pr.completed

        # TODO: awkward way to format time
        # removes microseconds from output
        return str(dt - self.created).split('.')[0]

    def is_manual(self):
        """Return true if this dataset is a manual input (rather than a file based input)"""
        return self.file_type == 'manual'

    def get_heatmap(self):
        """Return data in the heatmap format with embeded graphs"""
        output = {
            'rows': [],
            'cols': [],
        }
        # Track all locusts for all strains, Mutable variable populated by make_scatter!
        locusts = []
        for strain in self.strains.all():
            for _ in strain.get_prediction(locusts):
                pass # Populates locusts for square plots

        for strain in self.strains.all():
            row = {'name': strain.name, 'cols': []}
            output['rows'].append(row)

            for drug, (drprob, fpos, fneg, graph) in strain.get_prediction(locusts):
                if drug not in output['cols']:
                    output['cols'].append(drug)

                index = output['cols'].index(drug)
                row['cols'].extend([None] * (index - len(row['cols']) + 1))
                row['cols'][index] = {
                    'name': drug,
                    'scatter': graph,
                    'false_positive': fpos,
                    'false_negative': fneg,
                    'dr_probability': drprob,
                }

        return output

    def user_name(self):
        """Return the uploader's username"""
        if self.user:
            return self.user.username
        return 'n/a'

    def user_affiliation(self):
        """Return the uploader's addiliation"""
        if self.user:
            return self.user.affiliation
        return 'n/a'

    def user_email(self):
        """Return the uploader's email address"""
        if self.user:
            return self.user.email
        return 'n/a'

    def save(self, *args, **kwargs):
        """Override the save function to populate some fields"""
        if not self.id:
            super(PredictDataset, self).save(*args, **kwargs)

        if not self.md5:
            self.md5 = md5('%s%s' % (self.id, self.title)).hexdigest()

        if not self.file_directory:
            # We make a new directory in the BASE_DIR based on this object's
            # primary-key ID padded with zeros to make them fixed width.
            job_name = 'tbdata_' + str(self.id).zfill(8)
            self.file_directory = join(self.BASE_DIR, job_name)
            if not isdir(self.file_directory):
                os.makedirs(self.file_directory)

        return super(PredictDataset, self).save(*args, **kwargs)

    def get_full_json(self):
        """Returns the dataset as a serialised json"""
        return serializers.serialize('json', PredictDataset.objects.filter(id=self.id))


class PredictPipeline(Model):
    """Each file type can have possible pipelines, this provides a selection"""
    pipeline = ForeignKey(Pipeline, related_name='predicts')
    file_type = CharField(choices=PredictDataset.FILE_TYPES, max_length=25)
    is_default = BooleanField(default=False)

    def __str__(self):
        return "Pipeline %s" % str(self.pipeline)


class PredictStrain(Model):
    """Each strain uploaded for a dataset"""
    name = CharField(max_length=128)

    dataset = ForeignKey(PredictDataset, related_name='strains')
    pipeline = ForeignKey(Pipeline)
    piperun = ForeignKey(PipelineRun, null=True, blank=True, on_delete=SET_NULL)

    # We need two file slots for pair ended fastq files
    file_one = ForeignKey(UploadFile, null=True, blank=True,
                          related_name='link_a', on_delete=SET_NULL)
    file_two = ForeignKey(UploadFile, null=True, blank=True,
                          related_name='link_b', on_delete=SET_NULL)
    files = property(lambda self: [a for a in (self.file_one, self.file_two) if a])

    def run(self):
        """Runs this pipeline as set (even if run before)"""
        options = {'output_dir': self.dataset.file_directory}

        options['clean_files'] = []
        if self.file_one:
            options['file'] = self.file_one.fullpath
            options['file_one'] = self.file_one.fullpath
            if self.dataset.delete_sources:
                options['clean_files'].append(options['file_one'])
        if self.file_two:
            options['file'] = self.name
            options['file_two'] = self.file_two.fullpath
            if self.dataset.delete_sources:
                options['clean_files'].append(options['file_two'])

        name = slugify("{}.{}".format(self.dataset.title, self.name))
        self.piperun = self.pipeline.run(name, **options)
        self.save()
        return self.piperun.programs.filter(is_error=True).count() == 0

    @property
    def files(self):
        return [fh for fh in (self.file_one, self.file_two) if fh is not None]

    @property
    def status(self):
        """
        Returns the status of this strain prediction run. A combination of
        the files_status and run_status, with an extra status if a prediction
        file is detected or if the internal time limit has been reached.
        """
        files_status = self.files_status
        if files_status == STATUS_DONE:
            run_status = self.run_status
            if run_status == STATUS_DONE and self.has_prediction:
                return STATUS_READY
            if self.has_timedout():
                return STATUS_TIMEOUT
            return run_status + files_status
        if self.has_timedout():
            return STATUS_TIMEOUT
        return files_status

    def has_timedout(self):
        """Returns True if this prediction run has timed out"""
        return self.dataset.created < get_timeout()

    def get_status(self):
        """Return a string description of the status for this prediction strain"""
        return STATUS_LABELS[self.status]

    def get_status_level(self):
        """Return the bootstrap css color level for this status"""
        return STATUS_LEVELS[self.status]

    @property
    def files_status(self):
        """
        Returns
          - 0: Downloading not started yet
          - 1: The files are still downloading
          - 2: There was an error downloading the files.
          - 3: The files are ready for this pipeline
        """
        if not self.files:
            return STATUS_NOFILES
        for input_file in self.files:
            if input_file.retrieval_error:
                return STATUS_ERROR
            if not input_file.retrieval_start:
                return STATUS_WAIT
            if not input_file.retrieval_end:
                return STATUS_START
        return STATUS_DONE

    @property
    def run_status(self):
        """Return the pipeline run status without updating"""
        """
        Returns
          - 0: Pipeline not started yet
          - 1: Pipeline program is submitted or started
          - 2: Pipeline program is not complete yet
          - 3: The files are ready for this pipeline
        """
        qs = self.piperun.programs.all() if self.piperun else []
        if not qs:
            return STATUS_WAIT
        for program in qs:
            if program.is_error:
                return STATUS_ERROR
            if not program.is_submitted:
                return STATUS_WAIT
            if not program.is_complete:
                return STATUS_START
        return STATUS_DONE

    @property
    def has_prediction(self):
        """Returns true if prediction data is available"""
        return bool(self.prediction_file)

    @property
    def prediction_file(self):
        """Return a detected prediction file for this strain"""
        for (url, fn, name) in self.output_files:
            if fn.endswith('matrix.json'):
                return fn
        return None

    @property
    def has_lineage(self):
        """Returns true if lineage data is available"""
        return bool(self.lineage_file)

    @property
    def lineage_file(self):
        """Return a detected lineage file for this strain"""
        for (_, filename, _) in self.output_files:
            if filename.endswith('lineage.txt'):
                return filename
        return None

    @property
    def lineage(self):
        """The name of the lineage if possible"""
        filename = self.lineage_file
        if filename and os.path.isfile(filename):
            with open(filename, 'r') as fhl:
                data = fhl.read().strip('\n\r')
                if '\t' in data:
                    data = data.split('\t')
                    return {
                        'type': 'spoligo',
                        'spoligotype': data[1],
                        'match': data[-1],
                    }

                data = [lin.replace('lineage', '') for lin in data.split(',')]
                for x, lin in enumerate(data):
                    if x and not lin.startswith(data[x-1]):
                        return {
                            'type': 'mixed',
                            'all': data,
                        }

                return {
                    'type': 'lineage',
                    'first': data[0],
                    'lineage': data[-1],
                }
        return 'Not Found'

    @property
    def output_files(self):
        """Iterate over all output files and return media url and filename"""
        root_path = self.dataset.BASE_DIR
        root_url = join(settings.MEDIA_URL, 'data')
        if self.piperun:
            for run in self.piperun.programs.all():
                for fn in run.output_fn:
                    if os.path.isfile(fn):
                        url = None
                        if root_path in fn:
                            root_path = root_path.rstrip('/')
                            url = fn.replace(root_path, root_url)
                        yield (url, fn, basename(fn))

    def get_raw_prediction(self):
        """Get the raw data slightly bound better"""
        matrix_fn = self.prediction_file
        if matrix_fn and os.path.isfile(matrix_fn):
            with open(matrix_fn, 'r') as fhl:
                try:
                    (pr, m_A, m_B) = json.loads(fhl.read())
                except ValueError:
                    logging.error("Can't load prediction file: %s" % matrix_fn)
                    m_A = []

                # Rotate mutation matrix 90 degrees
                for name in m_A:
                    yield (name, zip(
                      [(b,c,d,e) for (a,b,c,d,e) in pr if a == name],
                      zip(*m_A[name]),
                      zip(*m_B[name]),
                    ))

    def get_prediction(self, locusts=None):
        """Get the prediction data formatted for heatmap and scatter plots"""
        for name, dat in self.get_raw_prediction():
            for (drug_code, dr, fp, fn), A, B in dat:
                try:
                    drug = Drug.objects.get(code__iexact=drug_code)
                except Drug.DoesNotExist:
                    sys.stderr.write("Can't find drug %s\n" % drug_code)
                    continue
                yield (drug_code, (dr, fp, fn, self.get_graph(drug, A, B, locusts)))

    def get_graph(self, drug, A, B, locusts=None):
        locusts = [] if locusts is None else locusts
        return [{
            "yAxis": "1",
            "cols": locusts,
            "key": key,
            "color": "rgba(%s)" % color,
            "values": list(self.make_scatter(locusts, datum)),
        } for key, color, datum in (
            ("Important", "255, 0, 0, 0.8", A),
            ("Other", "0, 0, 255, 0.17", B),
        )]

    def make_scatter(self, locusts, data):
        """Turn a mutations data in a by-gene plot"""
        regions = defaultdict(list)
        for gene in data:
            if gene:
                (_, region, _) = unpack_mutation_format(gene)
                regions[region].append(gene)
                if region not in locusts:
                    locusts.append(region)

        for x, locust in enumerate(locusts):
            ret = {"x": x, "y": 0, "size": 5, "tip": ["No mutations"]}
            if locust in regions:
                ret["y"] = len(regions[locust])
                ret["size"] = 9
                ret["tip"] = regions[locust]
            yield ret

    def __str__(self):
        return str(self.name)


class PredictDatasetNote(TimeStampedModel):
    """Notes of background processes"""
    dataset = ForeignKey(PredictDataset, related_name='notes')
    title = CharField(max_length=255)
    note = TextField()

    def __str__(self):
        return str(self.title)

    class Meta:
        ordering = ('-modified', '-created')
