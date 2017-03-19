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
Django forms for adding PredictDataset objects as well as a Confirmation form
"""
import os
import json

from django.utils.text import slugify
from django.core.urlresolvers import reverse_lazy
from django.forms import * 

from apps.dropbox.widgets import DropboxChooserWidget
from apps.dropbox.models import DropboxFile
from apps.mutations.fields import GeneticInputField
from apps.mutations.models import Mutation

from .models import PredictDataset, PredictStrain, PredictPipeline

FASTQ_FILES = ['.fastq', '.fastq.gz']
VCF_FILES = ['.vcf', '.vcf.gz']

class UploadForm(ModelForm):
    """
    Form for a user to enter a title, description, and dropbox files.
    """
    do_not_call_in_templates = True
    pipeline = ModelChoiceField(queryset=PredictPipeline.objects.all())
    my_status = PredictDataset.STATUS['DATASET_CONFIRMED']
    my_file_type = None
    ordered = 100

    class Meta:
        model = PredictDataset
        fields = ('title', 'description', 'status', 'user', 'file_type')
        widgets = {
          'status': HiddenInput(),
          'user': HiddenInput(),
          'file_type': HiddenInput(),
        }

    icon = classmethod(lambda cls: 'forms/%s.svg' % cls.my_file_type)
    enabled = classmethod(lambda cls: cls.pipeline_queryset().count() > 0)

    def __init__(self, *args, **kw):
        super(UploadForm, self).__init__(*args, **kw)
        pl = self.fields['pipeline']
        pl.queryset = self.pipeline_queryset()
        if pl.queryset.count() == 0:
            raise ValueError("No pipeline for type '%s'" % self.my_file_type)
        elif pl.queryset.count() == 1:
            pl.widget = HiddenInput()
            pl.initial = pl.queryset.get().pk
        else:
            try:
                pl.initial = pl.queryset.filter(is_default=True).get().pk
            except PredictPipeline.DoesNotExist:
                pass

    @classmethod
    def pipeline_queryset(cls):
        return PredictPipeline.objects.filter(file_type=cls.my_file_type)

    @classmethod
    def all_forms(cls):
        """Returns all the form classes inheriting this one"""
        subclasses = set()
        work = [cls]
        while work:
            parent = work.pop()
            for child in parent.__subclasses__():
                if child not in subclasses:
                    subclasses.add(child)
                    work.append(child)
        return sorted(subclasses, cmp=lambda a,b: cmp(a.ordered, b.ordered))

    def save(self, **kw):
        dataset = super(UploadForm, self).save(**kw)
        if not dataset.pk:
            return dataset
        for key, field in self.fields.items():
            if isinstance(field.widget, DropboxChooserWidget):
                for bucket_id, _, _ in field.widget.buckets:
                    self.save_dropbox(dataset, key + '_' + bucket_id)
                else:
                    self.save_dropbox(dataset, key)
        return dataset

    def get_strain_name(self, fn):
        """Gets the default name of the strain from the given filename"""
        if fn.endswith('.gz'):
            fn = fn[:-3]
        return ('file_one', fn.rsplit('.', 1)[0])

    def save_dropbox(self, dataset, key, files=None):
        if files is None:
            data = self.cleaned_data.get(key, self.data.get(key, '[]'))
            if data is None:
                return
            # Decode the data from the dropbox widget as json, this gives us
            # the file size, icon and filename which is very useful later.
            files = json.loads(data)

        for dropbox_file in files:
            (field, name) = self.get_strain_name(dropbox_file['name'])
            kw = {
              'pipeline': self.cleaned_data['pipeline'].pipeline,
            }
            (strain, _) = dataset.strains.get_or_create(name=name, defaults=kw)

            setattr(strain, field, DropboxFile.objects.create(
              name=key,
              file_directory=dataset.file_directory,
              filename=dropbox_file['name'],
              url=dropbox_file['link'],
              size=dropbox_file['bytes'],
              icon=dropbox_file['icon'],
            ))
            strain.save()


class ManualInputForm(UploadForm):
    my_file_type = PredictDataset.FILE_TYPE_MANUAL
    my_status = PredictDataset.STATUS['FILE_RETRIEVAL_SUCCESS']
    title = "Manually Entered Prediction"
    doc = "Create a prediction though manual selection of 1 or more mutations from a list. This option involves the shortest processing time but does assume that any non-entered mutations have been tested for and are absent."
    genetic_information = GeneticInputField(reverse_lazy('genes:json'))
    ordered = 20
    btn = 'default'

    def clean_genetic_information(self):
        data = self.cleaned_data.get('genetic_information')
        mutations = [m.strip() for m in data.split('\n') if m.strip()]
        name = slugify(self.cleaned_data.get('title'))
        (output, left_over) = Mutation.objects.matrix_csv(name, mutations)
        if left_over:
            raise ValidationError("Mutations not found in database: %s" % left_over)
        return output

    def save(self, *args, **kw):
        dataset = super(UploadForm, self).save(*args, **kw)
        if dataset and dataset.pk:
            data = self.cleaned_data.get('genetic_information')
            self.save_dropbox(dataset, 'manual', [
                {'name': 'matrix.csv', 'link': '-', 'bytes': len(data), 'icon': None}
            ])
            # We save the matrix file directly.
            dataset.strains.get().file_one.save_now(data)
        return dataset


class UploadVcfForm(UploadForm):
    my_file_type = PredictDataset.FILE_TYPE_VCF
    title = "Create VCF Prediction"
    doc = 'Minimal genotypic information for accurate resistance predictions are below. Genetic regions are listed in order of decreasing importance. For more detailed list of genetic variants see reference <a href="http://www.ncbi.nlm.nih.gov/pubmed/26910495">Farhat MR, Sultana R et al. Genetic Determinants of Drug Resistance in Mycobacterium tuberculosis and Their Diagnostic Value. AJRCCM 2016</a> and get <a href="https://en.wikipedia.org/wiki/Variant_Call_Format">more information</a> about the VCF format.'
    vcf_file = CharField(widget=DropboxChooserWidget(VCF_FILES), required=True,
        label="VCF Files", help_text="Variant Call Formated sequence data file. Multiple files can be selected, one vcf file per stain to compare.")
    ordered = 10
    btn = 'primary'


class UploadFastQSingleForm(UploadForm):
    my_file_type = PredictDataset.FILE_TYPE_FASTQ
    title = "Create FastQ Single-Ended Prediction"
    doc = "Create a prediction from a single-ended FastQ genetic sequence file. This option involves a large file and takes more time to process that the VCF or manual options."
    fastq_file = CharField(widget=DropboxChooserWidget(FASTQ_FILES), required=True,
        label="FastQ Files", help_text="FastQ files containing the single sequence read. Multiple files can be selected, one fastq file per strain to compare.")
    ordered = 5


class UploadFastQPairForm(UploadForm):
    my_file_type = PredictDataset.FILE_TYPE_FASTQ2
    title = "Create FastQ Pair-Ended Prediction"
    doc = "Create a prediction from a set of pair-ended FastQ genetic sequences. This option involves the largest files and takes more time to process that the VCF or manual options."
    fastq_file = CharField(widget=DropboxChooserWidget(FASTQ_FILES, buckets=[
        ('forward', "_R1.fastq _R1.fastq.gz", "Forward FastQ Files"),
        ('backward', "_R2.fastq _R2.fastq.gz", "Backward FastQ Files")]),
        required=True, label="FastQ Files",
        help_text="FastQ file containing the forward and backward sequence. Multiple strains can be selected for comparison.")
    ordered = 6

    def get_strain_name(self, fn):
        # We need to save fastq pair ended together, they come in here
        # in two buckets, so we save them in two fields (for now)
        name = fn.split('R1')[0].split('R2')[0].strip('_- ')
        return (['file_one', 'file_two']['R2' in fn], name)

    def clean_fastq_file(self):
        value = json.loads(self.cleaned_data['fastq_file'])
        if value:
            raise ValidationError("Unknown files were included, please select only accepted files.")
        r1 = set(self.clean_fastq_dir(self.data['fastq_file_forward'], 'R1'))
        r2 = set(self.clean_fastq_dir(self.data['fastq_file_backward'], 'R2'))
        extra = r1 ^ r2
        if extra:
            raise ValidationError("Unmatched files found: %s" % ", ".join(extra))

    def clean_fastq_dir(self, files, direction='R1'):
        for fastq_file in json.loads(files):
            name = fastq_file['name']
            key = "_%s.fastq" % direction
            if name.endswith('.gz'):
                name = name[:-3]
            if not name.endswith(key):
                raise ValidationError("Filename '%s' invalid, remember to include R1 or R2 suffix in each filename." % fastq_file['name'])
            yield name[:0-len(key)]

class NotificationForm(Form):
    success = BooleanField(required=False)
    result_data = CharField(widget=Textarea, required=False)

    def was_run_successful(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['success']

    def get_result_data(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['result_data']

