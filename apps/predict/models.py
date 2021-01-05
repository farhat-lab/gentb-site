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
import json
import logging

from hashlib import md5

import os
from os.path import join, isdir, basename

from collections import defaultdict
from datetime import timedelta

from model_utils.models import TimeStampedModel

from django.db.models import (
    Model, CASCADE, SET_NULL,
    ForeignKey, CharField, TextField, BooleanField, DecimalField, IntegerField, DateTimeField,
)
from django.conf import settings
from django.utils.text import slugify
from django.core import serializers
from django.urls import reverse
from django.utils.timezone import now
from django.utils.translation import ugettext_lazy as _

from apps.uploads.models import UploadFile
from apps.pipeline.models import Pipeline, PipelineRun, ProgramRun
from apps.mutations.models import Drug, GeneLocus

from .utils import lineage_spoligo, lineage_fast_caller, lineage_other_caller, filter_none

LOGGER = logging.getLogger('apps.predict')

# Basic status matrix, file_status + processing_status
STATUS = dict([
    ('FILE_NONE', (_('No Files Uploaded'), 'danger', 9)),
    ('FILE_WAIT', (_('Dataset Confirmed'), 'default', 10)),
    ('FILE_START', (_('File Retrieval Started'), 'primary', 8)),
    ('FILE_ERROR', (_('File Retrieval Failed'), 'danger', 10)),
    ('RUN_NONE', (_('No Strains to process'), 'danger', 10)),
    ('RUN_WAIT', (_('File Retrieval Success'), 'default', 8)),
    ('RUN_START', (_('Processing Started'), 'primary', 6)),
    ('RUN_ERROR', (_('Processing Failed'), 'danger', 4)),
    ('RUN_DONE', (_('Processing Success'), 'default', 10)),
    ('READY', (_('Prediction Ready'), 'success', 0)),
    ('INVALID', (_('Lacks Quality'),'warning', 8)),
    ('TIMEOUT', (_('Processing Timed Out'), 'danger', 3)),
])

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

    user = ForeignKey(settings.AUTH_USER_MODEL, related_name='datasets',
                      null=True, on_delete=SET_NULL)
    md5 = CharField(max_length=40, blank=True, db_index=True,\
            help_text='auto-filled on save')
    title = CharField('Dataset Title', max_length=255)
    file_type = CharField(choices=FILE_TYPES, max_length=25)
    delete_sources = BooleanField(default=False,\
        help_text="If this is checked, we will delete all your input files"\
            " downloaded from dropbox after running the predict.")
    has_notified = BooleanField(default=False,\
        help_text="Has this predict messaged the user about the status of their jobs.")

    description = TextField('Dataset description')
    file_directory = CharField(max_length=255, blank=True)

    class Meta:
        ordering = ('-created', 'title')

    def __str__(self):
        return str(self.title)

    def get_absolute_url(self):
        """Return a link to thedataset view"""
        if self.md5:
            return reverse('predict:view_single_dataset', kwargs=dict(slug=self.md5))
        return '/'

    def get_status(self):
        """Returns the status as a string"""
        return STATUS[self.status][0]

    def get_status_level(self):
        """Returns the btn/bootstrap color level for this status"""
        return STATUS[self.status][1]

    @property
    def status(self):
        """Returns the numeric level which this status has"""
        status = 'RUN_NONE'
        previous = 100
        for strain in self.strains.all():
            if STATUS[strain.status][2] < previous:
                status = strain.status
                previous = STATUS[status][2]
        return status

    @property
    def statuses(self):
        """Return a count of statuses"""
        ret = defaultdict(int)
        total = 0
        for strain in self.strains.all():
            ret[strain.status] += 1
            total += 1
        if not ret:
            ret['RUN_NONE'] = 1
        return [{
            'code': status,
            'count': count,
            'label': STATUS[status][0],
            'level': STATUS[status][1],
            'pc': "{:.2f}".format(count / total * 100),
        } for (status, count) in ret.items()]

    @property
    def directory_exists(self):
        """Returns true if the file_directory exists"""
        return os.path.isdir(self.file_directory)

    @property
    def has_prediction(self):
        """Returns true if any strain has a predict file"""
        return any(strain.has_prediction for strain in self.strains.all())

    @property
    def has_lineages(self):
        """Return true if any strain has a lineage file"""
        return any(strain.has_lineage for strain in self.strains.all())

    @property
    def has_output_files(self):
        """Return true if any strain has output files"""
        return any(list(strain.output_files) for strain in self.strains.all())

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
        for strain in self.strains.filter(results__isnull=True):
            strain.generate_results()

        strains = defaultdict(list)
        qset = PredictResult.objects.filter(strain__dataset=self, drug__isnull=False)
        vset = qset.values_list('id', 'strain__name', 'drug__code',\
                                'false_positive', 'false_negative', 'probability')

        for pk, strain, drug, fpos, fneg, prob in vset:
            cols = strains[strain]
            if drug not in output['cols']:
                output['cols'].append(drug)
                for ocol in strains.values():
                    if len(ocol) < len(output['cols']):
                        ocol.append(None)

            index = output['cols'].index(drug)
            cols.extend([None] * (index - len(cols) + 1))
            cols[index] = {
                'result_id': pk, 'name': drug,
                'false_positive': fpos, 'false_negative': fneg, 'dr_probability': prob,
            }
        for strain, cols in strains.items():
            output['rows'].append({'name': strain, 'cols': cols})
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
            key = '{}{}'.format(self.id, self.title)
            self.md5 = md5(key.encode('utf8')).hexdigest()

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

    def lineages(self):
        """Get a table of lineages"""
        header = []
        results = []
        for strain in self.strains.all():
            columns = [''] * len(header)
            for name, lineage in strain.lineages():
                if name in header:
                    columns[header.index(name)] = lineage
                else:
                    header.append(name)
                    columns.append(lineage)
            results.append({
                'strain': str(strain),
                'cols': columns,
            })
        for row in results:
            row['cols'] += [''] * (len(header) - len(row['cols']))
        return {
            'header': header,
            'rows': results,
        }


class PredictPipeline(Model):
    """Each file type can have possible pipelines, this provides a selection"""
    pipeline = ForeignKey(Pipeline, related_name='predicts', null=True, on_delete=SET_NULL)
    file_type = CharField(choices=PredictDataset.FILE_TYPES, max_length=25)
    is_default = BooleanField(default=False)

    def __str__(self):
        return "Pipeline %s" % str(self.pipeline)


class PredictStrain(Model):
    """Each strain uploaded for a dataset"""
    name = CharField(max_length=128)

    dataset = ForeignKey(PredictDataset, related_name='strains', on_delete=CASCADE)
    pipeline = ForeignKey(Pipeline, null=True, on_delete=SET_NULL)
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
        if files_status != 'FILE_DONE':
            return files_status

        run_status = self.run_status
        if run_status == 'RUN_DONE' and self.has_prediction:
            return 'READY'
        return run_status

    def has_timedout(self):
        """Returns True if this prediction run has timed out"""
        return self.dataset.created < get_timeout()

    def get_status(self):
        """Return a string description of the status for this prediction strain"""
        return STATUS[self.status][0]

    def get_status_level(self):
        """Return the bootstrap css color level for this status"""
        return STATUS[self.status][1]

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
            return 'FILE_NONE'
        for input_file in self.files:
            if input_file.retrieval_error:
                return 'FILE_ERROR'
            if not input_file.retrieval_start:
                return 'FILE_WAIT'
            if not input_file.retrieval_end:
                return 'FILE_START'
        return 'FILE_DONE'

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
            return 'RUN_WAIT'
        for program in qs:
            if program.job_state == 'TIMEOUT':
                return 'TIMEOUT'
            if program.is_error or program.job_state == 'INVALID':
                if program.program.quality_control:
                    return 'INVALID'
                return 'RUN_ERROR'
            if not program.is_submitted:
                return 'RUN_WAIT'
            if not program.is_complete:
                return 'RUN_START'
        return 'RUN_DONE'

    @property
    def has_prediction(self):
        """Returns true if prediction data is available"""
        return bool(self.prediction_file)

    @property
    def prediction_file(self):
        """Return a detected prediction file for this strain"""
        for (url, fn, name, life) in self.output_files:
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
        for (_, filename, _, life) in self.output_files:
            if filename.endswith('lineage.txt'):
                return filename
        return None

    def lineages(self):
        """The name of the lineage if possible"""
        filename = self.lineage_file
        if not filename:
            return []

        with open(filename, 'r') as fhl:
            data = fhl.read().strip('\n\r')

            # Very old format lineages
            if data and data[0] == '\t':
                return list(lineage_spoligo(data.split('\t')))
            if '\t' in data:
                return list(lineage_fast_caller(data))
            return list(lineage_other_caller(data))
        return []

    @property
    def output_files(self):
        """Iterate over all output files and return media url and filename"""
        root_path = self.dataset.BASE_DIR
        root_url = join(settings.MEDIA_URL, 'data')
        if self.piperun:
            for run in self.piperun.programs.all():
                life = run.output_life()
                for fn in run.output_fn:
                    if os.path.isfile(fn):
                        url = None
                        if root_path in fn:
                            root_path = root_path.rstrip('/')
                            url = fn.replace(root_path, root_url)
                        yield (url, fn, basename(fn), life)

    def get_raw_prediction(self):
        """Get the raw data slightly bound better"""
        matrix_fn = self.prediction_file
        if matrix_fn and os.path.isfile(matrix_fn):
            try:
                yield from self._prediction_from_file(matrix_fn)
            except Exception:
                yield None, False

    @staticmethod
    def _prediction_from_file(matrix_fn):
        m_A, m_B, m_C, m_D = {}, {}, {}, {}
        with open(matrix_fn, 'r') as fhl:
            parts = json.loads(fhl.read())
            (pr, m_A, m_B) = parts[:3]

            # this is different because it assumes that there's only one strain.
            name = list(m_A)[0]
            m_C[name] = [([None] * len(m_A[name][0]))] * len(m_A[name])
            m_D[name] = [([None] * len(m_A[name][0]))] * len(m_A[name])

            ex_1, ex_2 = parts[-2:] if len(parts) == 5 else ({}, {})
            for index, value in ex_1.items():
                m_C[name][int(index)] = filter_none(value)
            for index, value in ex_2.items():
                m_D[name][int(index)] = filter_none(value)

        # Rotate mutation matrix 90 degrees
        for name in m_A:
            yield (name, zip(
                [(b,c,d,e) for (a,b,c,d,e) in pr if a == name],
                zip(*m_A[name]),
                zip(*m_B[name]),
                zip(*m_C[name]),
                zip(*m_D[name]),
            ))

    def generate_results(self):
        """Populate the database from the matrix files"""
        self.results.all().delete()

        for _, dat in self.get_raw_prediction():
            if dat is False: # Error in raw prediction getting, prevent asking again
                self.results.create(drug=None)
                break

            for (drug_code, dr, fneg, fpos), *data in dat:
                try:
                    drug = Drug.objects.get(code__iexact=drug_code)
                except Drug.DoesNotExist:
                    continue
                res, _ = self.results.get_or_create(drug=drug, defaults={
                    'probability':dr, 'false_positive':fpos, 'false_negative':fneg})

                for cat, datum in enumerate(data):
                    for mutation in datum:
                        if not mutation:
                            continue
                        try:
                            locus = GeneLocus.objects.for_mutation_name(mutation, True)
                        except ValueError:
                            continue
                        (obj, created) = res.loci.get_or_create(category=cat + 1, locus=locus,
                                                                 defaults={'mutations': mutation})
                        if not created:
                            obj.mutations += "\n" + mutation

    def __str__(self):
        return str(self.name)


class PredictResult(Model):
    """A resulting prediction of a specific drug and strain"""
    strain = ForeignKey(PredictStrain, related_name='results', on_delete=CASCADE)
    drug = ForeignKey(Drug, related_name='strain_predictions', null=True, blank=True, on_delete=CASCADE)

    updated = DateTimeField(auto_now=True)

    false_negative = DecimalField("False Negative Rate", null=True, blank=True, decimal_places=5, max_digits=10)
    false_positive = DecimalField("False Postive Rate", null=True, blank=True, decimal_places=5, max_digits=10)
    probability = DecimalField("Drug Resistance Probability", null=True, blank=True, decimal_places=5, max_digits=10)

    class Meta:
        unique_together = [('strain', 'drug')]

class PredictResultLocus(Model):
    """A specific locus for the result"""
    result = ForeignKey(PredictResult, related_name='loci', on_delete=CASCADE)
    locus = ForeignKey(GeneLocus, related_name='prediction_results', on_delete=CASCADE)
    category = IntegerField(choices=(
        (0, 'Unknown'),
        (1, 'Important'),
        (2, 'Other'),
        (3, 'New'),
        (4, 'Lineage SNPs'),
    ))
    mutations = TextField(default='')

    class Meta:
        unique_together = ['result', 'locus', 'category']

class PredictDatasetNote(TimeStampedModel):
    """Notes of background processes"""
    dataset = ForeignKey(PredictDataset, related_name='notes', on_delete=CASCADE)
    title = CharField(max_length=255)
    note = TextField()

    def __str__(self):
        return str(self.title)

    class Meta:
        ordering = ('-modified', '-created')
