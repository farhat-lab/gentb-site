import os
import json

import logging
LOGGER = logging.getLogger('apps.predict')

from hashlib import md5
from os.path import join, isfile, isdir, basename

from collections import defaultdict

from django.db.models import *
from model_utils.models import TimeStampedModel

from django.conf import settings
from django.utils.text import slugify
from django.core import serializers
from django.core.urlresolvers import reverse
from django.utils.timezone import now
from django.utils.translation import ugettext_lazy as _

from apps.dropbox.models import DropboxFile
from apps.pipeline.models import Pipeline, PipelineRun
from apps.mutations.models import Drug
from apps.mutations.utils import unpack_mutation_format

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

    user = ForeignKey(settings.AUTH_USER_MODEL, related_name='datasets')
    md5 = CharField(max_length=40, blank=True, db_index=True,
            help_text='auto-filled on save')
    title = CharField('Dataset Title', max_length=255)
    file_type = CharField(choices=FILE_TYPES, max_length=25)

    description = TextField('Dataset description')
    file_directory = CharField(max_length=255, blank=True)

    status = PositiveIntegerField(default=1, choices=STATUS_CHOICES)
    is_error = property(lambda self: self.status in [0, 4, 8])
    is_busy = property(lambda self: self.status in [3, 6])

    def __str__(self):
        return str(self.title)

    def get_absolute_url(self):
        return reverse('predict:view_single_dataset', kwargs=dict(slug=self.md5))

    @property
    def media_url(self):
        """Return the location of the file directory accessable via url"""
        return join(settings.MEDIA_URL, 'data', basename(self.file_directory.rstrip('/')))

    @property
    def directory_exists(self):
        """Returns true if the file_directory exists"""
        return os.path.isdir(self.file_directory)

    @property
    def has_prediction(self):
        return isfile(join(self.file_directory, 'output', 'matrix.json'))

    def is_manual(self):
        return self.file_type == 'manual'

    def make_scatter(self, locusts, data):
        regions = defaultdict(list)
        for gene in data:
            if gene:
                (index, region, mutation) = unpack_mutation_format(gene)
                regions[region].append(gene)


        for x, locust in enumerate(locusts):
            ret = {"x": x, "y": 0, "size": 5, "tip": ["No mutations"]}
            if locust in regions:
                ret["y"] = len(regions[locust])
                ret["size"] = 9
                ret["tip"] = regions[locust]
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
                    
                    locusts = list(drug.mutations.values_list('gene_locus__name', flat=True).distinct())
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

    def save(self, *args, **kwargs):
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


class PredictPipeline(Model):
    """Each file type can have possible pipelines, this provides a selection"""
    pipeline = ForeignKey(Pipeline)
    file_type = CharField(choices=PredictDataset.FILE_TYPES, max_length=25)
    is_default = BooleanField(default=False)

    def __str__(self):
        return "Pipeline %s" % str(self.pipeline)


class PredictStrain(Model):
    """Each strain uploaded for a dataset"""
    name = CharField(max_length=128)

    dataset = ForeignKey(PredictDataset, related_name='strains')
    pipeline = ForeignKey(Pipeline)
    piperun = ForeignKey(PipelineRun, null=True, blank=True)
    
    # We need two file slots for pair ended fastq files
    file_one = ForeignKey(DropboxFile, null=True, blank=True, related_name='link_a')
    file_two = ForeignKey(DropboxFile, null=True, blank=True, related_name='link_b')
    files = property(lambda self: [a for a in (self.file_one, self.file_two) if a])

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

