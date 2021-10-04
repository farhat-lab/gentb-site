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

from django.utils.text import slugify
from django.urls import reverse_lazy
from django.core.exceptions import ValidationError
from django.forms import (
    ModelForm, Form, ModelChoiceField, HiddenInput, Textarea,
    BooleanField, CharField,
)

from apps.uploads.fields import UploadField
from apps.uploads.models import UploadFile
from apps.mutations.fields import GeneticInputField
from apps.mutations.models import Mutation

from .models import PredictDataset, PredictStrain, PredictPipeline

FASTQ_FILES = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']
VCF_FILES = ['.vcf', '.vcf.gz']

class UploadForm(ModelForm):
    """
    Form for a user to enter a title, description, and dropbox files.
    """
    do_not_call_in_templates = True
    pipeline = ModelChoiceField(queryset=PredictPipeline.objects.filter(is_retired=False))
    my_file_type = None
    ordered = 100

    class Meta:
        model = PredictDataset
        fields = ('title', 'description', 'user', 'file_type', 'delete_sources')
        widgets = {
          'user': HiddenInput(),
          'file_type': HiddenInput(),
        }

    icon = classmethod(lambda cls: 'forms/%s.svg' % cls.my_file_type)
    enabled = classmethod(lambda cls: cls.pipeline_queryset(is_default=True).count() > 0)

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
        self.uploads = [key for key, field in self.fields.items()
            if isinstance(field, UploadField)]

    @classmethod
    def pipeline_queryset(cls, **kw):
        return PredictPipeline.objects.filter(file_type=cls.my_file_type, is_retired=False, **kw)

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
        return sorted(subclasses, key=lambda a: a.ordered)

    @classmethod
    def get_form(cls, frm_id):
        for item in cls.all_forms():
            if item.my_file_type == frm_id:
                return item
        raise KeyError("Form not found: %s" % frm_id)

    def save(self, **kw):
        dataset = super(UploadForm, self).save(**kw)
        if dataset.pk:
            for key in self.uploads:
                self.save_upload(dataset, self.cleaned_data.get(key, {}))
        return dataset

    def clean(self):
        """Clean the form"""
        data = super(UploadForm, self).clean()
        strain_fields = [f.name for f in PredictStrain._meta.get_fields()]
        for key in self.uploads:
            for fl in data.get(key, {}):
                if fl is not None and fl not in strain_fields:
                    raise ValidationError("Upload '%s' is invalid." % fl)
        return data

    def save_upload(self, dataset, field):
        for bucket, upload_files in field.items():
            for upload_file in upload_files:
                (strain, _) = dataset.strains.get_or_create(
                    name=upload_file.name,
                    defaults={
                      'pipeline': self.cleaned_data['pipeline'].pipeline,
                    }
                )
                upload_file.conclude_upload(dataset.file_directory, dataset.user)
                attr = bucket or 'file_one'
                if not hasattr(strain, attr):
                    raise AttributeError("File Bucket '%s' doesn't exist!" % attr)
                setattr(strain, attr, upload_file)
                strain.save()


class ManualInputForm(UploadForm):
    my_file_type = PredictDataset.FILE_TYPE_MANUAL
    doc_title = "Manually Entered Prediction"
    doc_template = "predict/types/manual.html"
    doc = "Create a prediction though the manual selection of 1 or more mutations from a list."
    genetic_information = GeneticInputField(reverse_lazy('genes:json'))
    ordered = 20
    btn = 'default'

    def __init__(self, *args, **kw):
        super(ManualInputForm, self).__init__(*args, **kw)
        self.fields['delete_sources'].widget = HiddenInput()

    def clean_genetic_information(self):
        data = self.cleaned_data.get('genetic_information')
        mutations = [m.strip() for m in data.split('\n') if m.strip()]
        name = slugify(self.cleaned_data.get('title'))
        (output, left_over) = Mutation.objects.matrix_csv(name, mutations)
        if left_over:
            raise ValidationError("Mutations are not available in prediction: %s" % ", ".join(list(left_over)))
        return output

    def save(self, *args, **kw):
        """
        Manual save by-passes any upload file in favour of constructing
        it's own matrix csv file and adding it directly to the strain.
        """
        dataset = super(UploadForm, self).save(*args, **kw)
        if dataset and dataset.pk:
            data = self.cleaned_data.get('genetic_information')

            # Create an UploadFile without the Input Field
            matrix = UploadFile.objects.create(
                name='matrix',
                file_directory=dataset.file_directory,
                filename='matrix.csv',
                size=len(data))

            # We save the matrix file directly
            matrix.save_now(data)
            self.save_upload(dataset, {'file_one': [matrix]})
        return dataset


class UploadVcfForm(UploadForm):
    my_file_type = PredictDataset.FILE_TYPE_VCF
    doc_title = "Varient Call File (VCF)"
    doc_template = "predict/types/vcf.html"
    doc = "Create a prediction from a variant call file in VCF format."
    vcf_file = UploadField(extensions=VCF_FILES, required=True,
        label="VCF Files", help_text="Variant Call Formated sequence data file. Multiple files can be selected, one vcf file per stain to compare.")
    ordered = 10
    btn = 'primary'


class UploadFastQSingleForm(UploadForm):
    my_file_type = PredictDataset.FILE_TYPE_FASTQ
    doc_title = "Single-Ended FastQ File"
    doc_template = "predict/types/fastq_single.html"
    doc = "Create a prediction from single-ended FastQ genetic sequence."
    fastq_file = UploadField(extensions=FASTQ_FILES, required=True,
        label="FastQ Files", help_text="FastQ files containing the single sequence read. Multiple files can be selected, one fastq file per strain to compare.")
    ordered = 5


class UploadFastQPairForm(UploadForm):
    my_file_type = PredictDataset.FILE_TYPE_FASTQ2
    doc_title = "Pair-Ended FastQ File"
    doc_template = "predict/types/fastq_pair.html"
    doc = "Create a prediction from a set of pair-ended FastQ genetic sequences."
    fastq_file = UploadField(extensions=FASTQ_FILES, buckets=[
        ('file_one', "^(.+)[\._\- ][Rr]?1\.(fastq|fq)(?:\.gz)?$", "Forward FastQ Files", 'file_two'),
        ('file_two', "^(.+)[\._\- ][Rr]?2\.(fastq|fq)(?:\.gz)?$", "Backward FastQ Files", 'file_one')],
        required=True, label="FastQ Files",
        help_text="FastQ file containing the forward and backward sequence. Multiple strains can be selected for comparison.")
    ordered = 6


class NotificationForm(Form):
    success = BooleanField(required=False)
    result_data = CharField(widget=Textarea, required=False)

    def was_run_successful(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['success']

    def get_result_data(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['result_data']

