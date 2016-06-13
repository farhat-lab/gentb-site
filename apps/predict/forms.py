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
                url = self.cleaned_data.get(key, None)
                if url:
                    dataset.files.create(name=key, url=url)
        return dataset


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
    vcf_file = CharField(widget=DropboxChooserWidget(['.vcf']), required=True,
        label="VCF File", help_text="Variant Call Formated sequence data file")


class UploadFastQPairForm(UploadForm):
    fastq_one = CharField(widget=DropboxChooserWidget(['.fastq']), required=True,
        label="Forward Read", help_text="FastQ file containing the forward read sequence.")
    fastq_two = CharField(widget=DropboxChooserWidget(['.fastq']), required=True,
        label="Reverse Read", help_text="FastQ file containing the reverse read sequence.")


class UploadFastQSingleForm(UploadForm):
    fastq_one = CharField(widget=DropboxChooserWidget(['.fastq']), required=True,
        label="FastQ File", help_text="FastQ file containing the single sequence read.")


class NotificationForm(Form):
    success = BooleanField(required=False)
    result_data = CharField(widget=Textarea, required=False)

    def was_run_successful(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['success']

    def get_result_data(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['result_data']

