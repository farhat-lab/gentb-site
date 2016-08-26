"""
Django forms for adding PredictDataset objects as well as a Confirmation form

"""
import os
import json

from django.utils.text import slugify
from django.core.urlresolvers import reverse_lazy
from django.forms import * 

from apps.dropbox.widgets import DropboxChooserWidget
from apps.mutations.fields import GeneticInputField
from apps.mutations.models import Mutation

from .models import PredictDataset

FASTQ_FILES = ['.fastq', '.fastq.gz']
VCF_FILES = ['.vcf', '.vcf.gz']

class UploadForm(ModelForm):
    """
    Form for a user to enter a title, description, and dropbox files.
    """
    class Meta:
        model = PredictDataset
        fields = ('title', 'description', 'status', 'user', 'file_type', 'fastq_type')
        widgets = {
          'status': HiddenInput(),
          'user': HiddenInput(),
          'file_type': HiddenInput(),
          'fastq_type': HiddenInput(),
        }

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

    def save_dropbox(self, dataset, key):
        data = self.cleaned_data.get(key, self.data.get(key, '[]'))
        if data is None:
            return
        files = json.loads(data)
        for dropbox_file in files:
            dataset.files.create(
              name=key,
              filename=dropbox_file['name'],
              url=dropbox_file['link'],
              size=dropbox_file['bytes'],
              icon=dropbox_file['icon'],
            )


class ManualInputForm(UploadForm):
    """
    Manually enter genetic information for prediction.
    """
    genetic_information = GeneticInputField(reverse_lazy('genes:json'))

    def clean_genetic_information(self):
        data = self.cleaned_data.get('genetic_information')
        mutations = [m.strip() for m in data.split('\n') if m.strip()]
        name = slugify(self.cleaned_data.get('title'))
        (output, left_over) = Mutation.objects.matrix_csv(name, mutations)
        if left_over:
            raise ValidationError("Mutations not found in database: %s" % left_over)
        return output

    def save(self, *args, **kw):
        obj = super(UploadForm, self).save(*args, **kw)
        if obj and obj.pk:
            path = os.path.join(obj.file_directory, 'output')
            if not os.path.isdir(path):
                os.makedirs(path)
            with open(os.path.join(path, 'matrix.csv'), 'w') as fhl:
                fhl.write(self.cleaned_data.get('genetic_information'))
        return obj


class UploadVcfForm(UploadForm):
    vcf_file = CharField(widget=DropboxChooserWidget(VCF_FILES), required=True,
        label="VCF Files", help_text="Variant Call Formated sequence data file. Multiple files can be selected, one vcf file per stain to compare.")


class UploadFastQPairForm(UploadForm):
    fastq_file = CharField(widget=DropboxChooserWidget(FASTQ_FILES, buckets=[
        ('forward', "_R1.fastq _R1.fastq.gz", "Forward FastQ Files"),
        ('backward', "_R2.fastq _R2.fastq.gz", "Backward FastQ Files")]),
        required=True, label="FastQ Files",
        help_text="FastQ file containing the forward and backward sequence. Multiple strains can be selected for comparison.")

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


class UploadFastQSingleForm(UploadForm):
    fastq_file = CharField(widget=DropboxChooserWidget(FASTQ_FILES), required=True,
        label="FastQ Files", help_text="FastQ files containing the single sequence read. Multiple files can be selected, one fastq file per strain to compare.")


class NotificationForm(Form):
    success = BooleanField(required=False)
    result_data = CharField(widget=Textarea, required=False)

    def was_run_successful(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['success']

    def get_result_data(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['result_data']

